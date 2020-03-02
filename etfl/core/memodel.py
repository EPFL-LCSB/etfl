# -*- coding:utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

Core for the ME-part

"""

import numpy as np
import optlang
import pandas as pd
import sympy
from cobra import Model, Reaction, Gene
from cobra.core import Solution, DictList
from collections import defaultdict, OrderedDict
from Bio.SeqUtils import molecular_weight

from tqdm import tqdm

from ..utils.parsing import parse_gpr
from ..utils.utils import replace_by_enzymatic_reaction, replace_by_me_gene
from .genes import ExpressedGene
from .dna import DNA
from .rna  import mRNA,rRNA, tRNA
from .enzyme import Enzyme, Peptide
from .reactions import EnzymaticReaction, ProteinComplexation, \
    TranslationReaction, TranscriptionReaction, DegradationReaction, DNAFormation
from .expression import build_trna_charging, enzymes_to_gpr_no_stoichiometry, \
    make_stoich_from_aa_sequence, make_stoich_from_nt_sequence, \
    degrade_peptide, degrade_mrna, _extract_trna_from_reaction
from ..optim.constraints import CatalyticConstraint, ForwardCatalyticConstraint,\
    BackwardCatalyticConstraint, EnzymeMassBalance, \
    rRNAMassBalance, mRNAMassBalance, tRNAMassBalance, DNAMassBalance, \
    GrowthCoupling, TotalCapacity, ExpressionCoupling, EnzymeRatio, \
    GrowthChoice, EnzymeDegradation, mRNADegradation, \
    SynthesisConstraint, SOS1Constraint,\
    InterpolationConstraint, RNAPAllocation
from ..optim.variables import ModelVariable, GrowthActivation, \
    EnzymeVariable, RibosomeUsage, RNAPUsage, \
    FreeEnzyme, BinaryActivator, InterpolationVariable, DNAVariable, \
    GrowthRate, GenericVariable

from .allocation import add_dummy_expression,add_dummy_mrna,add_dummy_peptide,\
    add_dummy_protein, add_interpolation_variables, \
    MRNA_WEIGHT_CONS_ID, PROT_WEIGHT_CONS_ID, \
    MRNA_WEIGHT_VAR_ID, PROT_WEIGHT_VAR_ID, \
    DNA_WEIGHT_CONS_ID, DNA_WEIGHT_VAR_ID, DNA_FORMATION_RXN_ID, \
    define_prot_weight_constraint, define_mrna_weight_constraint, \
    define_dna_weight_constraint, get_dna_synthesis_mets

from pytfa.core.model import LCSBModel
from pytfa.optim.reformulation import linearize_product
from pytfa.optim import LinearizationConstraint
from pytfa.optim.utils import chunk_sum, symbol_sum
from pytfa.utils.logger import get_bistream_logger
from pytfa.utils.str import camel2underscores
from pytfa.optim.utils import copy_solver_configuration


class MEModel(LCSBModel, Model):
    def __init__(self, model=Model(), name = None,
                 growth_reaction='',
                 mu_range = None, n_mu_bins = 1,
                 big_M = 1000,
                 *args, **kwargs):

        """


        :param model: The input model
        :type model: cobra.Model
        :param mu:  (Facultative) Mean growth rate to constraint the model
        :param mu_error: (Facultative) Absolute error on mu to constraint the model
        :type mu_error: float > 0
        :param mu_range: (Facultative) Min-Max growth rate to constraint the model
        :type mu_range: tuple (l,u)
        :param n_mu_bins: (Facultative) In how many intervals to separate the
                        growth rate for the linearization
        :param args:
        :param kwargs:

        """

        name = 'ETFL_' + name if name is not None else 'ETFL_model'

        LCSBModel.__init__(self, model, name)

        self.logger = get_bistream_logger('ME model' + str(self.name))
        self.parent = model
        if model is not None:
            self.sanitize_varnames()

        self.init_etfl(big_M, growth_reaction, mu_range,
                       n_mu_bins, name)

    def init_etfl(self, big_M, growth_reaction, mu_range, n_mu_bins, name):

        self.big_M = big_M
        self._var_dict = dict()
        self._cons_dict = dict()
        self.logger.info('# ETFL Model {} initialized'.format(name))
        self._growth_reaction_id = growth_reaction
        self._mu_range = mu_range
        self._n_mu_bins = n_mu_bins
        if mu_range is not None:
            self._mu = self.add_variable(kind=GrowthRate,
                                         hook=self,
                                         id_='total',  # Will read MU_total
                                         lb=mu_range[0],
                                         ub=mu_range[1])
            self.init_mu_variables()
        else:
            # message = """ You need to supply mu_range."""

            message = "Empty model initialized"
            # raise ValueError(message)
            self.logger.info(message)
        self.aa_dict = dict()
        self.rna_nucleotides = dict()
        self.trna_dict = dict()
        self.dilution_terms = dict()
        self.enzymes = DictList()
        self.mrnas = DictList()
        self.rrnas = DictList()
        self.trnas = DictList()
        self.peptides = DictList()
        self.transcription_reactions = DictList()
        self.translation_reactions = DictList()
        self.complexation_reactions = DictList()
        self.degradation_reactions = DictList()

        self.dna = None
        self.ribosome = OrderedDict()
        self.rnap = OrderedDict()

        self.coupling_dict = dict()

    @property
    def mu(self):
        return self._mu

    @property
    def mu_max(self):
        return self._mu_range[-1]

    # @mu.setter
    # def mu(self, val, epsilon = None):
    #     if epsilon is None:
    #         epsilon = self.solver.configuration.tolerances.feasibility
    #
    #     self._mu.lb = val-epsilon
    #     self._mu.ub = val+epsilon

    def make_mu_bins(self):
        from numpy import linspace
        bounds = linspace(self.mu.variable.lb, self.mu.variable.ub, self._n_mu_bins+1)
        bins = zip(bounds[:-1], bounds[1:])
        self.mu_bins = tuple(((x[0]+x[1])/2, x) for x in bins)


    @property
    def n_mu_bins(self):
        return len(self.mu_bins)

    def init_mu_variables(self):
        """
        Necessary for the zeroth order approximation of mu:

        .. math::

            mu \in [0.1, 0.9] , nbins = 8
            => mu = 0.15 OR mu = 0.25 OR ... OR mu = 0.85

        Using binary expansion of the bins instead of a list of 0-1s
        described `here <https://orinanobworld.blogspot.ch/2013/07/integer-variables-and-quadratic-terms.html>`_

        :return:
        """

        self.make_mu_bins()
        ga = list()
        N = self.n_mu_bins
        n_vars = np.int(np.ceil(np.log2(N)))

        for e in range(n_vars):
            ga.append(self.add_variable(kind=GrowthActivation,
                                        hook=self,
                                        id_=str(2 ** e)))

        # Force that only one growth range can be chosen:
        # b0*2^0 + b1*2^1 + b2*2^2 + ... + bn*2^n <= n_bins

        choice_expr = sum(ga)
        self.add_constraint(kind=GrowthChoice,
                            hook=self,
                            expr=choice_expr,
                            id_='growth',
                            ub=self.n_mu_bins,
                            lb=0)

        # Couple growth
        v_fwd = self.growth_reaction.forward_variable
        v_bwd = self.growth_reaction.reverse_variable

        # |v_net - mu| <= bin_width
        bin_half_width = max([(x[1] - x[0]) / 2 for _, x in self.mu_bins])

        the_integer = symbol_sum([(2 ** i) * ga_i for i, ga_i in enumerate(ga)])

        binarized_mu = self.mu.variable.lb + the_integer * self.mu_approx_resolution

        growth_coupling_expr = v_fwd - v_bwd - binarized_mu

        self.add_constraint(kind=GrowthCoupling,
                            hook=self.growth_reaction,
                            expr=growth_coupling_expr,
                            ub=bin_half_width,
                            lb=-1 * bin_half_width)


        # So that the solver spends less time looking for an ub on the objective
        # when optimizing for growth
        self.growth_reaction.upper_bound = self.mu.variable.ub + self.mu_approx_resolution

        # Update the variable indices
        self.regenerate_variables()
        self.regenerate_constraints()

    @property
    def mu_approx_resolution(self):
        return (self.mu.variable.ub - self.mu.variable.lb) / self.n_mu_bins

    @property
    def growth_reaction(self):
        """
        Returns the growth reaction of the model. Useful because tied to the
        growth variable
        
        :return:
        """
        if self._growth_reaction_id:
            return self.reactions.get_by_id(self._growth_reaction_id)
        else:
            return None

    @growth_reaction.setter
    def growth_reaction(self, reaction_id):
        """
        The growth_reaction is set by supplying the id of the candidate reaction

        :param reaction_id: an id within the model
        :type reaction_id: str
        :return:
        """
        rxn = self.reactions.get_by_id(reaction_id)
        self._growth_reaction_id = rxn.id

    def add_nucleotide_sequences(self, sequences):
        """

        :param sequences:
        :return:
        """

        for gene_id, seq in sequences.items():
            if gene_id in self.genes:
                new = replace_by_me_gene(self, gene_id, seq)

            else:
                self.logger.warning('Model has no gene {}, Adding it'.format(gene_id))
                new = ExpressedGene(id= gene_id, name = gene_id, sequence=seq)
                self.add_genes([new])

            self._make_peptide_from_gene(gene_id)

    def _make_peptide_from_gene(self, gene_id):
        free_pep = Peptide(id=gene_id,
                           name='Peptide, {}'.format(gene_id),
                           gene_id=gene_id)
        free_pep._model = self
        self.peptides += [free_pep]

    def add_dummies(self, nt_ratios, mrna_kdeg, mrna_length, aa_ratios,
                    enzyme_kdeg, peptide_length):
        """

        Create dummy peptide and mrna to enforce mrna and peptide production.
        Can be used to account for the missing data for all mrnas and proteins.

        :param nt_ratios:
        :param mrna_kdeg:
        :param mrna_length:
        :param aa_ratios:
        :param enzyme_kdeg:
        :param peptide_length:
        :param gtp:
        :param gdp:
        :param h2o:
        :param h:
        :return:
        """

        add_interpolation_variables(self)

        # Create a dummy gene and override the sequences with input data
        dummy_sequence = 'N'*mrna_length
        dummy_gene = ExpressedGene(id='dummy_gene',
                                   name='Dummy Gene',
                                   sequence=dummy_sequence)
        dummy_gene._rna = dummy_sequence
        dummy_gene._peptide = 'X'*peptide_length

        self.add_genes([dummy_gene])

        dummy_mrna = add_dummy_mrna(self, dummy_gene, mrna_kdeg, mrna_length, nt_ratios)

        dummy_peptide = add_dummy_peptide(self, aa_ratios, dummy_gene, peptide_length)

        dummy_protein = add_dummy_protein(self, dummy_peptide, enzyme_kdeg)

        add_dummy_expression(self, aa_ratios, dummy_gene, dummy_peptide, dummy_protein,
                                   peptide_length)

    def add_essentials(self, essentials, aa_dict, rna_nucleotides,
                       rna_nucleotides_mp):
        """
        Marks important metabolites for expression

        :param essentials: A dictionary of important metabolites to met id

            **Example :**

            .. code-block:: python

                essentials = {
                            'atp': 'atp_c',
                            'adp': 'adp_c',
                            'amp': 'amp_c',
                            ...
                            'h2o': 'h2o_c',
                            'h': 'h_c'}
                        }

        :param aa_dict: A dictionary of aminoacid letter to amicoacid met id

            **Example :**

            .. code-block:: python

                aa_dict = {
                            'A':'ala__L_c',
                            'R':'arg__L_c',
                            ...
                        }

        :param rna_nucleotides: A dictionary of RNA nucleotide triphosphate
            letter to nucleotideTP met id

            **Example :**

            .. code-block:: python

                rna_nucleotides = {
                            'A':'atp_c',
                            'U':'utp_c',
                            ...
                        }

        :param rna_nucleotides_mp: A dictionary of RNA nucleotide monophosphate
            letter to nucleotideMP met id

            **Example :**

            .. code-block:: python

                rna_nucleotides_mp = {
                            'A':'amp_c',
                            'U':'ump_c',
                            ...
                        }

        :return:
        """

        self.essentials = essentials
        self.aa_dict = aa_dict
        self.rna_nucleotides = rna_nucleotides
        self.rna_nucleotides_mp = rna_nucleotides_mp


    def build_expression(self):
        """
        Given a dictionary from amino acids nucleotides to metabolite names,
        goes through the list of genes in the model that have sequence
        information to build transcription and translation reactions

        :return:
        """

        aa_dict = self.aa_dict

        atp = self.essentials['atp']
        amp = self.essentials['amp']
        ppi = self.essentials['ppi']
        h2o = self.essentials['h2o']
        h = self.essentials['h']

        self.trna_dict = build_trna_charging(self,aa_dict,atp,amp,ppi,h2o,h)
        self.add_trnas([item for sublist in self.trna_dict.values()
                        for item in sublist if isinstance(item, tRNA)])

        # Check that the ribosomes have been added
        if not self.ribosome:
            raise Exception(
                'A ribosome has to be added with the add_ribosome method')

        # Check that the RNAP has been added
        if not self.rnap:
            raise Exception(
                'A RNA Polymerase has to be added with the add_rnap method')

        self.express_genes(self.genes)

    def express_genes(self, gene_list):
        """
        Adds translation and transcription reaction to the genes in the
        provided list

        :param gene_list:
        :type gene_list: Iterable of str or ExpressedGene
        :return:
        """

        for gene in gene_list:

            if isinstance(gene, str):
                the_gene = self.genes.get_by_id(gene)
            else:
                the_gene = self.genes.get_by_id(gene.id)

            if not isinstance(the_gene, ExpressedGene):
                continue

            # Build the translation
            self._add_gene_translation_reaction(the_gene)
            # Build the transcription
            self._add_gene_transcription_reaction(the_gene)

    def _add_gene_translation_reaction(self, gene):
        """

        :param gene: A gene of the model that has sequence data
        :type gene: ExpressedGene
        :return:
        """

        gtp = self.essentials['gtp']
        gdp = self.essentials['gdp']
        pi  = self.essentials['pi']
        h2o = self.essentials['h2o']
        h   = self.essentials['h']

        if gene.translated_by is None:
            translating_enzymes = self.ribosome.values()
        else:
            translating_enzymes = [self.enzymes.get_by_id(x)
                                   for x in gene.translated_by]

        rxn = TranslationReaction(
            id='{}_translation'.format(gene.id),
            name='Translation, {}'.format(gene.id),
            gene_id= gene.id,
            enzymes=translating_enzymes,
            upper_bound=1,
            scaled=True)

        self.add_reactions([rxn])

        aa_stoichiometry = make_stoich_from_aa_sequence(gene.peptide,
                                                        self.aa_dict,
                                                        self.trna_dict,
                                                        gtp,
                                                        gdp,
                                                        pi,
                                                        h2o,
                                                        h
                                                        )

        _extract_trna_from_reaction(aa_stoichiometry, rxn)

        rxn.add_metabolites(aa_stoichiometry, rescale=True)

        free_peptide = self.peptides.get_by_id(gene.id)

        rxn.add_peptide(free_peptide)

        # Add ribosome as necessary enzyme
        rxn.gene_reaction_rule = enzymes_to_gpr_no_stoichiometry(rxn)
        self.translation_reactions += [rxn]

    def _add_gene_transcription_reaction(self, gene):
        """
        Adds the transcription reaction related to a gene

        :param gene: A gene of the model that has sequence data
        :type gene: ExpressedGene
        :return:
        """

        ppi  = self.essentials['ppi']

        nt_stoichiometry = make_stoich_from_nt_sequence(gene.rna,
                                                        self.rna_nucleotides,
                                                        ppi
                                                        )

        if gene.transcribed_by is None:
            transcribing_enzymes = self.rnap.values()
        else:
            transcribing_enzymes = [self.enzymes.get_by_id(x)
                                   for x in gene.transcribed_by]

        rxn = TranscriptionReaction(
            id=self._get_transcription_name(gene.id),
            name='Transcription, {}'.format(gene.id),
            gene_id= gene.id,
            enzymes=transcribing_enzymes,
            upper_bound=1,
            scaled=True)
        self.add_reactions([rxn])

        rxn.add_metabolites(nt_stoichiometry)

        # Add rnap as necessary enzyme
        rxn.gene_reaction_rule = enzymes_to_gpr_no_stoichiometry(rxn)


        self.transcription_reactions += [rxn]

    def add_trna_mass_balances(self):
        """
        Once the tRNAs, transcription and translation reactions have been added,
        we need to add the constraints:

        d/dt [charged_tRNA]   =  v_charging - sum(nu_trans*v_trans) - mu*[charged_tRNA]
        d/dt [uncharged_tRNA] = -v_charging + sum(nu_trans*v_trans) - mu*[uncharged_tRNA]

        The stoichiometries are set from the reaction dict in _extract_trna_from_reaction

        We also need to scale the tRNAs in mRNA space and unscale the translation:

        .. code::

            d/dt σ_m * [*charged_tRNA] =    +- σ_m * v_charging
                                            -+ σ_m/σ_p*sum(nu_tsl*σ_p*v_tr)
                                            -  σ_m * mu*[*charged_tRNA]

            d/dt [*charged_tRNA]_hat =      +- σ_m * v_charging
                                            -+ σ_m/σ_p * sum( nu_tsl * v_tr_hat)
                                        -  mu*[*charged_tRNA]_hat

        :return:
        """

        translation_fluxes = self.translation_reactions.list_attr('net')

        for _, (charged_trna, uncharged_trna, charging_rxn) in self.trna_dict.items():

            # Charged tRNAs are generated with the charging reaction, consumed
            # by translation
            charged_stoichs = [translation.trna_stoich[charged_trna.id] for
                               translation in self.translation_reactions]

            # ys are negative because charged tRNAs are consumed in translation
            v_tsl_c = symbol_sum(
                [x * y
                 for x, y in zip(charged_stoichs, translation_fluxes)])

            charged_expr = charging_rxn.forward_variable \
                           + v_tsl_c
            self.add_mass_balance_constraint(synthesis_flux=charged_expr,
                                             macromolecule=charged_trna,
                                             queue=True)

            # Uncharged tRNAs are generated whenever translation happens,
            # consumed by charging
            uncharged_stoichs = [translation.trna_stoich[uncharged_trna.id] for
                                 translation in self.translation_reactions]

            v_tsl_u = symbol_sum(
                [x * y
                 for x, y in zip(uncharged_stoichs, translation_fluxes)])

            uncharged_expr = -1 * charging_rxn.forward_variable \
                             + v_tsl_u

            self.add_mass_balance_constraint(synthesis_flux=uncharged_expr,
                                             macromolecule=uncharged_trna,
                                             queue=True)

        self._push_queue()

    def add_enzymatic_coupling(self, coupling_dict):
        """
        Couples the enzymatic reactions maximal rates with the Enzyme
        availability
        The coupling dictionary looks like:

        .. code-block:: python

            coupling_dict : {
                            'reaction_id_1':[   enzyme_instance_1,
                                                enzyme_instance_2],
                            'reaction_id_2':[   enzyme_instance_3,
                                                enzyme_instance_4,
                                                enzyme_instance_5],

        :param coupling_dict: A dictionary of reaction ids to enzyme lists
        :type coupling_dict: {str:list(Enzyme)}
        :return:
        """
        self.coupling_dict.update(coupling_dict)
        self.add_enzymes(coupling_dict.values())

        # /!\ We modify the reaction list
        # self.add_gene_reactions()

        # Generic reactions <-> Enzymes coupling
        for rid in tqdm(coupling_dict, desc='cat. constraints'):
            r = self.reactions.get_by_id(rid)

            if isinstance(r, EnzymaticReaction) and r.id in coupling_dict:
                # This is a proper enzymatic reaction and we can directly apply
                # the constraint
                self.logger.debug('Applying catalytic constraint to {}'. \
                                  format(rid))
                r.add_enzymes(coupling_dict[r.id])
                self.apply_enzyme_catalytic_constraint(r)
            elif not isinstance(r, EnzymaticReaction) and r.id in coupling_dict:
                # This reaction needs to be transformed to an EnzymaticReaction
                self.logger.debug('Transforming and applying catalytic constraint to {}'. \
                                  format(rid))
                #TODO : Add enzymatic_reaction dictlist ??
                enzymes = coupling_dict[r.id]
                enz_r = replace_by_enzymatic_reaction(self, r.id, enzymes, scaled=False)
                self.apply_enzyme_catalytic_constraint(enz_r)
            else:
                self.logger.error('Could not find reaction {} in the coupling dictionary'.format(r.id))

        # update variable and constraints attributes
        self._push_queue()
        self.regenerate_constraints()
        self.regenerate_variables()

    def apply_enzyme_catalytic_constraint(self, reaction):
        """
        Apply a catalytic constraint using a gene-enzymes reaction rule (GPR)

        :param reaction:
        :return:
        """

        v_max_fwd = dict()
        v_max_bwd = dict()

        # Write v_max constraint
        fwd_variable = reaction.forward_variable
        bwd_variable = reaction.reverse_variable

        for e, enz in enumerate(reaction.enzymes):
            # If the enzymes has the same kcat for both directions
            # v_fwd  <= kcat_fwd [E]
            # v_fwd - kcat_fwd [E] <= 0

            v_max_fwd[e] = enz.kcat_fwd * enz.concentration
            v_max_bwd[e] = enz.kcat_bwd * enz.concentration

        # Formulating the scaling factor on the max kcat
        k_f = max([x.kcat_fwd for x in reaction.enzymes])
        k_b = max([x.kcat_bwd for x in reaction.enzymes])
        # k_f = np.median([x.kcat_fwd for x in self.enzymes])
        # k_b = np.median([x.kcat_bwd for x in self.enzymes])
        E_m = max([x.scaling_factor for x in reaction.enzymes])

        # v_fwd <= sum(kcat_i*E_i)
        # for all i, E_i <= E_max (= 1g/gDW)
        # v_fwd / sum(kcat_i*E_i^max) <= sum(kcat_i*E_i) / sum(kcat_i*E_i^max) (<= 1)

        enz_constraint_expr_fwd = (fwd_variable - sum(v_max_fwd.values()))/(k_f*E_m)
        enz_constraint_expr_bwd = (bwd_variable - sum(v_max_bwd.values()))/(k_b*E_m)


        self.add_constraint(kind=ForwardCatalyticConstraint, hook=reaction,
                            expr=enz_constraint_expr_fwd, ub=0, queue=True)
        self.add_constraint(kind=BackwardCatalyticConstraint, hook=reaction,
                            expr=enz_constraint_expr_bwd, ub=0, queue=True)

    def add_mass_balance_constraint(self, synthesis_flux, macromolecule=None,
                                    queue=False):
        """
        Adds a mass balance constraint of the type

        .. math::

            d[E]/dt = 0 <=> v_synthesis - k_deg*[M] - μ*[M] = 0

        for a macromolecule (mRNA or enzyme)

        :param synthesis_flux:
        :param macromolecule:
        :return:
        """

        kwargs = dict()

        try:
            z = self.dilution_terms[macromolecule.id]
            self.logger.warning('Dilution term already exists for {}'.format(macromolecule.id))
        except KeyError:
            z = self.linearize_me(macromolecule, queue=queue)

        if isinstance(macromolecule, Enzyme):
            kind = EnzymeMassBalance
            hook = macromolecule
            # expr = synthesis_flux.scaled_net \
            #        - macromolecule.degradation.scaled_net \
            #        - 1/macromolecule.kdeg * z
                   # - 1/self.mu_max * z
            expr = synthesis_flux.net \
                   - macromolecule.degradation.net \
                   - macromolecule.scaling_factor * z
            # expr = 1/macromolecule.degradation.scaling_factor * \
            #        ( synthesis_flux.net \
            #        - macromolecule.degradation.net \
            #        - macromolecule.scaling_factor * z)
        elif isinstance(macromolecule, mRNA):
            kind = mRNAMassBalance
            hook = macromolecule
            # sigma = 1/macromolecule.degradation.scaling_factor
            # expr = sigma * synthesis_flux.net \
            #        - macromolecule.degradation.scaled_net \
            #        - 1/macromolecule.kdeg * z
                   # - 1/self.mu_max * z
            expr = synthesis_flux.net \
                   - macromolecule.degradation.net \
                   - macromolecule.scaling_factor * z

            # expr = 1/macromolecule.degradation.scaling_factor * \
            #        ( synthesis_flux.net \
            #        - macromolecule.degradation.net \
            #        - macromolecule.scaling_factor * z)
        elif isinstance(macromolecule, tRNA):
            kind = tRNAMassBalance
            kwargs['id_'] = macromolecule.id
            hook = self
            # mu_ub = self.mu_max
            # sigma = mu_ub * macromolecule.scaling_factor
            # The synthesis flux comes unscaled [flux units]
            # expr = 1/sigma * synthesis_flux - 1/mu_ub * z
            expr = synthesis_flux - macromolecule.scaling_factor * z
        elif isinstance(macromolecule, DNA):
            kind = DNAMassBalance
            kwargs['id_'] = 'dna'
            hook = self
            # mu_ub = self.mu_max
            # expr = synthesis_flux.scaled_net - 1/mu_ub *z
            expr = synthesis_flux.net - macromolecule.scaling_factor * z
        else:
            raise Exception('Macromolecule type not recognized: {}'
                            .format(macromolecule))

        self.add_constraint(kind=kind,
                            hook=hook,
                            expr=expr,
                            lb=0, ub=0,
                            queue = queue,
                            **kwargs)

    def linearize_me(self, macromolecule, queue=False):
        """
        Performs Petersen linearization on μ*E to keep a MILP problem

        :return:
        """

        E = macromolecule.variable

        # ga_i is a binary variable for the binary expansion f the fraction on N
        # of the max growth rate
        ga_vars = self.get_ordered_ga_vars()

        out_expr = self.mu.variable.lb

        # Build z =   ga_0*2^0*mu_max/N * [E]
        #           + ga_1*2^1*mu_max/N * [E]
        #           + ...
        #           + ga_n*2^n*mu_max/N * [E]

        for i, ga_i in enumerate(ga_vars):
            # Linearization step for ga_i * [E]
            model_z_i = linearize_product(model = self, b = ga_i, x=E, queue=queue)

            out_expr += (2 ** i) * self.mu_approx_resolution * model_z_i

        self.dilution_terms[macromolecule.id] = out_expr

        return out_expr

    def get_ordered_ga_vars(self):
        """
        Returns in order the variables that discretize growth
        :return:
        """
        # ga_i is a binary variable for the binary expansion f the fraction on N
        # of the max growth rate
        ga_vars = self.get_variables_of_type(GrowthActivation)
        ga_vars = sorted(ga_vars, key=lambda x: x.ix)
        return ga_vars


    def _prep_enzyme_variables(self, enzyme):
        """
        Reads Enzyme.composition to find complexation reaction from enzyme information

        :param reaction:
        :type reaction: cobra.Reaction
        :return:
        """

        #1. Complexation

        # Does this enzyme already have a complexation reaction ?
        # This happens if an enzyme is used in several reactions
        if enzyme.complexation is not None:
            complexation = enzyme.complexation
        else:
            complexation = self.make_enzyme_complexation(enzyme)

        enzyme.init_variable()

        #2. Also add degradation
        self._add_enzyme_degradation(enzyme, scaled=True, queue=True)

        #3. Finally make the mass balance
        # Cannot queue, if the same enzyme is to be added twice
        self.add_mass_balance_constraint(complexation, enzyme, queue=False)

    def make_enzyme_complexation(self, enzyme):
        """
        Makes the complexation reaction and attached it to its enzyme

        :param enzyme:
        :return:
        """
        if not enzyme.composition:
            self.logger.warning('Enzyme {} has no composition'
                                .format(enzyme.id))
            return None

        this_id = '{}_complex'.format(enzyme.id)
        this_name = '{} Complexation'.format(enzyme.id)
        complexation = ProteinComplexation(id=this_id,
                                                name=this_name,
                                                target=enzyme,
                                                # upper_bound=1,
                                                scaled=True)

        try:
            peptides = {self.peptides.get_by_id(k): -v \
                        for k, v in enzyme.composition.items()}
        except KeyError:
            missing_genes = '.'.join(enzyme.composition.keys())
            self.logger.warning('No nucleotide sequence found for '
                                'some of these genes {}'.format(missing_genes))
            return None

        self.add_reactions([complexation])

        complexation.add_peptides(peptides)

        # Post processing
        self.complexation_reactions+= [complexation]
        enzyme.complexation = complexation

        return complexation

    def add_enzymes(self, enzyme_list, prep = True):
        """
        Adds an Enzyme object, or iterable of Enzyme objects, to the model
        :param enzyme_list:
        :type enzyme_list:Iterable(Enzyme) or Enzyme
        :param prep: whether or not to add complexation, degradation, and mass
            balance constraints (needs to be overridden for dummies for example)
        :type prep: Boolean
        :return:
        """
        if not hasattr(enzyme_list, '__iter__'):
            enzyme_list = [enzyme_list]
        else:
            enzyme_list = list(enzyme_list)
        if len(enzyme_list) == 0:
            return None

        # unpacking
        if not isinstance(enzyme_list[0],Enzyme):
            enzyme_list = [x for item in enzyme_list for x in item]


        # First check whether the enzymes exist in the model
        # Also enzymes could be declared twice for different reactions,
        # hence turn the list into a set
        enzyme_list = [x for x in set(enzyme_list) if x.id not in self.enzymes]

        enz_ids = [x.id for x in enzyme_list]

        tot_ids = len(enz_ids)
        unique_ids = len(set(enz_ids))
        if tot_ids != unique_ids:
            msg = '{} duplicate enzyme IDs detected'.format(tot_ids-unique_ids)
            self.logger.error(msg)
            raise KeyError(msg)

        for enz in enzyme_list:
            enz._model = self

        self.enzymes += enzyme_list

        if prep:
            for enz in tqdm(enzyme_list, desc='enz. vars'):
                self._prep_enzyme_variables(enz)


    def add_mrnas(self, mrna_list, add_degradation=True):
        """
        Adds a mRNA object, or iterable of mRNA objects, to the model
        :param mrna_list:
        :type mrna_list:Iterable(mRNA) or mRNA
        :return:
        """
        if not hasattr(mrna_list, '__iter__'):
            mrna_list = [mrna_list]
        if len(mrna_list) == 0:
            return None

            # First check whether the mRNAs exist in the model
            mrna_list = [x for x in mrna_list if x.id not in self.mrnas]

        for mrna in mrna_list:
            mrna._model = self
            mrna.init_variable()

            if add_degradation:
                self._add_mrna_degradation(mrna, scaled=True, queue=True)

        self.mrnas += mrna_list
        self._push_queue

    def add_trnas(self, trna_list):
        """
        Adds a tRNA object, or iterable of tRNA objects, to the model
        :param trna_list:
        :type trna_list:Iterable(tRNA) or tRNA
        :return:
        """
        if not hasattr(trna_list, '__iter__'):
            trna_list = [trna_list]
        if len(trna_list) == 0:
            return None

        # First check whether the tRNAs exist in the model
        trna_list = [x for x in trna_list if x.id not in self.trnas]

        for trna in trna_list:
            trna._model = self
            trna.init_variable()

        self.trnas += trna_list

    def add_dna(self, dna):
        """
        Adds a DNA object to the model

        :param dna:
        :type dna: DNA
        :return:
        """

        dna._model = self
        dna.init_variable()

        self.dna = dna

    def remove_enzymes(self, enzyme_list):
        """
        Removes an Enzyme object, or iterable of Enzyme objects, from the model

        :param enzyme_list:
        :type enzyme_list:Iterable(Enzyme) or Enzyme
        :return:
        """
        if not hasattr(enzyme_list, '__iter__'):
            enzyme_list = [enzyme_list]
        if len(enzyme_list) == 0:
            return None

        # First check whether the metabolites exist in the model
        enzyme_list = [x for x in enzyme_list if x.id not in self.enzymes]

        for enz in enzyme_list:
            self.remove_reactions(enz.degradation)
            self.remove_reactions(enz.complexation)
            self.enzymes.pop(enz.id)


    def _add_enzyme_degradation(self, enzyme, scaled=True, queue=False):
        """
        Given an enzyme, adds the corresponding degradation reaction

        :param enzyme:
        :type enzyme: Enzyme
        :param scaled: Indicates whether scaling should be performed (see manuscript)
        :type scaled: bool
        :param queue: Indicates whether to add the variable directly or
                        in the next batch
        :type queue: bool
        :return:
        """

        h2o = self.essentials['h2o']

        if enzyme.kdeg is None or np.isnan(enzyme.kdeg):
            return None

        complex_dict = enzyme.complexation.metabolites
        deg_stoich = defaultdict(int)
        for peptide, stoich in complex_dict.items():
            the_pep = self.peptides.get_by_id(peptide.id)
            degradation_mets = degrade_peptide(the_pep,
                                               self.aa_dict,
                                               h2o)
            for k,v in degradation_mets.items():
                deg_stoich[k]+=-1*v*stoich # stoich is negative

        self._make_degradation_reaction(deg_stoich,
                                        enzyme,
                                        EnzymeDegradation,
                                        scaled=scaled,
                                        queue=queue)


    def _add_mrna_degradation(self, mrna, scaled = True, queue=False):
        """
        Given an mRNA, adds the corresponding degradation reaction

        :param mrna:
        :type mrna: mRNA
        :param scaled: Indicates whether scaling should be performed (see manuscript)
        :type scaled: bool
        :param queue: Indicates whether to add the variable directly or
                        in the next batch
        :type queue: bool
        :return:
        """

        h2o = self.essentials['h2o']
        h = self.essentials['h']

        if mrna.kdeg is None or np.isnan(mrna.kdeg):
            return None

        degradation_mets = degrade_mrna(mrna, self.rna_nucleotides_mp, h2o, h)

        self._make_degradation_reaction(degradation_mets,mrna,mRNADegradation,
                                        scaled=scaled, queue=queue)


    def _make_degradation_reaction(self, deg_stoich, macromolecule,
                                   kind, scaled, queue=False):
        """
        given a degradation stoichiometry, makes the corresponding degradation
        reaction

        :param deg_stoich: stoichiometry of the degradation
        :type deg_stoich: dict({:class:`cobra.core.Species:Number`})
        :param macromolecule: the macromalecule being degraded. Used for binding
                                the degradation constraint
        :type macromolecule: Macromolecule
        :param kind: kind of constraint
        :type kind: mRNADegradation or EnzymeDegradation
        :param scaled: Indicates whether scaling should be performed (see manuscript)
        :type scaled: bool
        :param queue: Indicates whether to add the variable directly or
                        in the next batch
        :type queue: bool
        :return:
        """

        reaction = DegradationReaction(id='{}_degradation'.format(macromolecule.id),
                                       macromolecule=macromolecule,
                                       scaled=scaled)

        if scaled:
            reaction.upper_bound = 1

        # Assignment to model must be done before since met dict kas string keys
        self.add_reactions([reaction])
        self.degradation_reactions += [reaction]
        reaction.add_metabolites(deg_stoich, rescale = True)
        # Couple with the expression constraint v_deg = k_deg [E]
        # Scaled into v_deg_hat = E_hat
        # expr = reaction.scaled_net \
        #        - (macromolecule.kdeg / self.mu_max) * macromolecule.scaled_concentration
        # expr = reaction.scaled_net - macromolecule.scaled_concentration
        expr = reaction.net - macromolecule.kdeg * macromolecule.concentration
        self.add_constraint(kind=kind,
                            hook=macromolecule,
                            expr=expr,
                            lb=0,
                            ub=0,
                            queue=queue)


    def populate_expression(self, constrain_polymerases=True):
        """
        Populates RNAP, ribosomes, and their number on their respective templates

        :return:
        """
        self._populate_rnap()
        self._populate_ribosomes()
        self._push_queue()

        for the_mrna in tqdm(self.mrnas, desc='populating expression'):
            self._constrain_polysome(the_mrna)

        if self.dna is not None and constrain_polymerases:
            for the_gene in tqdm(self.genes, desc='constraining transcription'):
                self._constrain_polymerase(the_gene)

        self._push_queue()
        self.regenerate_variables()
        self.regenerate_constraints()

    def _constrain_polysome(self, the_mrna):
        """
        Add the coupling between mRNA availability and ribosome charging
        The number of ribosomes assigned to a mRNA species is lower than
        the number of such mRNA times the max number of ribosomes that can sit
        on the mRNA:
        [RPi] <= loadmax_i*[mRNAi]

        loadmax is : len(peptide_chain)/size(ribo)
        "Their distance from one another along the mRNA is at least the size
        of the physical footprint of a ribosome (≈20 nm, BNID 102320, 100121)
        which is the length of about 60 base pairs (length of
        nucleotide ≈0.3 nm, BNID 103777), equivalent to ≈20 aa."
        "http://book.bionumbers.org/how-many-proteins-are-made-per-mrna-molecule/"

        Hence:
        [RPi] <= L_nt/Ribo_footprint * [mRNA]

        :return:
        """
        ribo_footprint_size = 60  # [bp] see docstring
        rib_usage_vars = self.get_variables_of_type(RibosomeUsage)
        # Get the synthesis_flux
        syn_id = self._get_transcription_name(the_mrna.id)
        syn = self.transcription_reactions.get_by_id(syn_id)

        # TODO: Move this to another function
        # Add the mass balance constraint for the mrna
        self.add_mass_balance_constraint(syn, the_mrna, queue=True)

        # Get the ribosomes assigned to this translation
        RPi_hat = rib_usage_vars.get_by_id(the_mrna.id)

        # Get the mRNA concentration
        mrna_hat = the_mrna.scaled_concentration

        polysome_size = len(the_mrna.gene.rna) / ribo_footprint_size

        # σ_m is the mRNA scaling factor,
        # σ_r is the ribosome scaling factor
        # [RPi] <= Lmrna/Lrib * [mRNAi]
        # [RPi]_hat  <= Lmrna/Lrib * σ_m/σ_r * [mRNAi]_hat
        #

        # nondimensionalized:
        scaling_factor = the_mrna.scaling_factor / RPi_hat.scaling_factor
        expression_coupling = RPi_hat \
                              - polysome_size * scaling_factor * mrna_hat
        # expression_coupling = RPi_hat - polysome_size * mrna_hat

        # Add expression coupling
        self.add_constraint(kind=ExpressionCoupling,
                            hook=the_mrna.gene,
                            expr=expression_coupling,
                            queue=True,
                            ub=0)

    def _constrain_polymerase(self, the_gene):
        """
        Add the coupling between DNA availability and RNAP charging
        The number of RNAP assigned to a gene locus is lower than
        the number of such loci times the max number of RNAP that can sit
        on the locus:
        [RNAPi] <= loadmax_i*[# of loci]*[DNA]

        loadmax is : len(nucleotide chain)/size(RNAP)

        "The footprint of RNAP II [...] covers approximately 40 nt and is
        nearly symmetrical [...]."
        BNID    107873
        Range 	~40 Nucleotides

        Hence:
        [RNAPi] <= loadmax_i*[# of loci]*[DNA]

        :return:
        """

        rnap_usage_vars = self.get_variables_of_type(RNAPUsage)

        if not the_gene in rnap_usage_vars:
            self.logger.debug('Gene {} has no RNAPUsage var. '
                              'No RNAP allocation constraint added.'
                              .format(the_gene.id))
            return None

        rnap_alloc = self._get_rnap_allocation_expression(the_gene)

        # Add expression coupling
        self.add_constraint(kind=RNAPAllocation,
                            hook=the_gene,
                            expr=rnap_alloc,
                            queue=True,
                                ub=0)

    def _get_rnap_allocation_expression(self, the_gene):
        """
        Given a gene and the appropritate RNAP Usage variable, returns the
        expression of the RNAP allocation constraint

        :param the_gene:
        :param RNAPi_hat:
        :return:
        """
        rnap_usage_vars = self.get_variables_of_type(RNAPUsage)

        # Get the RNAP assigned to this transcription
        RNAPi_hat = rnap_usage_vars.get_by_id(the_gene.id)

        # Get the number of loci
        n_loci = the_gene.copy_number

        rnap_footprint_size = 40  # [bp] see docstring of ```_constrain_polymerases'''
        loadmax = len(the_gene.sequence) / rnap_footprint_size

        # σ_dna is the DNA scaling factor,
        # σ_rp is the RNAP scaling factor
        # [RNAPi]      <= Lgene/Lrnap              * n_loci * [DNA]
        # [RNAPi]_hat  <= Lgene/Lrnap * σ_dna/σ_rp * n_loci * [DNA]_hat
        #
        # Non-dimensionalized

        scaling_factor = self.dna.scaling_factor / RNAPi_hat.scaling_factor

        rnap_alloc = RNAPi_hat \
                     - loadmax * n_loci * scaling_factor * self.dna.scaled_concentration

        return rnap_alloc

    def edit_gene_copy_number(self, gene_id):
        """
        Edits the RNAP allocation constraints if the copy number of a gene
        changes.

        :param gene_id:
        :return:
        """
        the_gene = self.genes.get_by_id(gene_id)

        rnap_alloc = self.get_constraints_of_type(RNAPAllocation)
        rnap_usage = self.get_variables_of_type(RNAPUsage)
        if gene_id in rnap_alloc and gene_id in rnap_usage:
            # We need to edit the variable
            # Modify the polymerase constraint
            cons = rnap_alloc.get_by_id(gene_id)
            new_expr = self._get_rnap_allocation_expression(the_gene)
            cons.change_expr(new_expr)
            self.logger.debug('Changed RNAP allocation for gene {}'.format(gene_id))
        else:
            # There is no variable, maybe the allocation constraints
            # have not been made yet.
            # self.model._add_rnap_allocation(rnap_usage, self)
            pass

    def recompute_translation(self):
        """

        :return:
        """

        for rib_id in self.ribosome:
            #1.1 Find the RNAP capacity constraint
            cons = self.get_constraints_of_type(TotalCapacity).get_by_id(rib_id)

            #1.2 Edit the RNAP capacity constraint
            new_total_capacity = self._get_rib_total_capacity()
            cons.change_expr(new_total_capacity)

    def recompute_transcription(self):
        """

        :return:
        """

        # TODO: No need to recompute for the one we just added.
        # Keep tabs ? use an except list
        for rnap_id in self.rnap:
            # 1.1 Find the RNAP capacity constraint
            cons = self.get_constraints_of_type(TotalCapacity).get_by_id(rnap_id)

            # 1.2 Edit the RNAP capacity constraint
            new_total_capacity = self._get_rnap_total_capacity()
            cons.change_expr(new_total_capacity)

    def recompute_allocation(self):
        """

        :return:
        """

        #0. Does the model have allocation constraints ?
        interpolation_constraints = self.get_constraints_of_type(InterpolationConstraint)
        interpolation_variables = self.get_variables_of_type(InterpolationVariable)
        if not interpolation_constraints:
            return None


        #1. Remove previous allocation constraints

        # mRNA
        mrna_weight_def_cons = interpolation_constraints.get_by_id(MRNA_WEIGHT_CONS_ID)
        mrna_weight_var = interpolation_variables.get_by_id(MRNA_WEIGHT_VAR_ID)

        # Proteins
        prot_weight_def_cons = interpolation_constraints.get_by_id(PROT_WEIGHT_CONS_ID)
        prot_weight_var = interpolation_variables.get_by_id(PROT_WEIGHT_VAR_ID)

        #DNA
        dna_weight_def_cons  = interpolation_constraints.get_by_id(DNA_WEIGHT_CONS_ID)
        dna_weight_var = interpolation_variables.get_by_id(DNA_WEIGHT_VAR_ID)

        for the_cons in [mrna_weight_def_cons, prot_weight_def_cons, dna_weight_def_cons]:
            self.remove_constraint(the_cons)

        #2. Apply new allocation constraints
        define_prot_weight_constraint(self,prot_weight_var)
        define_mrna_weight_constraint(self,mrna_weight_var)

        chromosome_len = self.dna.len
        gc_content= self.dna.gc_ratio
        define_dna_weight_constraint(self, self.dna, dna_weight_var, gc_content, chromosome_len)

    def _get_transcription_name(self, the_mrna_id):
        """
        Given an mrna_id, gives the id of the corresponding transcription reaction
        :param the_mrna_id:
        :type the_mrna_id: str
        :return: str
        """
        return '{}_transcription'.format(the_mrna_id)

    def _get_translation_name(self, the_peptide_id):
        """
        Given an mrna_id, gives the id of the corresponding translation reaction
        :param the_peptide_id:
        :type the_peptide_id: str
        :return: str
        """
        return '{}_translation'.format(the_peptide_id)

    def get_translation(self, the_peptide_id):
        """
        Given an peptide_id, gives the translation reaction
        :param the_peptide_id:
        :type the_peptide_id: str
        :return: TranslationReaction
        """
        return self.reactions.get_by_id(self._get_translation_name(the_peptide_id))

    def get_transcription(self, the_peptide_id):
        """
        Given an mrna_id, gives corresponding transcription reaction
        :param the_mrna_id:
        :type the_mrna_id: str
        :return: TranscriptionReaction
        """
        return self.reactions.get_by_id(self._get_transcription_name(the_peptide_id))

    def add_rnap(self, rnap, free_ratio=0):
        """
        Adds the RNA Polymerase used by the model.

        :param rnap:
        :type rnap: Ribosome
        :return:
        """

        if rnap.id in self.rnap:
            raise KeyError('A RNAP with this id ({}) alread exists in '
                           'the model'.format(rnap.id))
        else:
            self.rnap[rnap.id] = rnap

        self.add_enzymes(rnap)

        if free_ratio > 0:
            self._add_free_enzyme_ratio(rnap, free_ratio)


    def _populate_rnap(self):
        """
        Once RNAP have been assigned to the model, we still need to link
        them to the rest of the variables and constraints. This function creates
        the mass balance constraint on the RNAP, as well as the total
        RNAP capacity constraint
        :return:
        """

        # v_complexation =   complexation.forward_variable  \
        #                  - complexation.reverse_variable

        # Create the mass balance constraint
        # 1 -> Write the RNAP mass balance
        # for this_rnap in self.rnap.values():
        #     self.add_mass_balance_constraint(this_rnap.complexation, this_rnap)
        # UPDATE: We now do this in add_rnap :)

        # 2 -> Parametrize all the transcription reactions with RNAP vmax
        for trans_rxn in self.transcription_reactions:
            self.apply_rnap_catalytic_constraint(trans_rxn, queue=True)
        self._push_queue()

        # 3 -> Add RNAP capacity constraint

        self.regenerate_variables()

        usage = self._get_rnap_total_capacity()

        # usage = (sum_RMs + self.RNAPf[this_rnap.id].unscaled - this_rnap.concentration) / \
        #         this_rnap.scaling_factor

        # Create the capacity constraint
        self.add_constraint(kind=TotalCapacity,
                            hook=self,
                            id_ = 'rnap',
                            expr=usage,
                            lb = 0,
                            ub = 0,
                            )

        # update variable and constraints attributes
        self.regenerate_constraints()
        self.regenerate_variables()

    def _get_rnap_total_capacity(self):
        all_rnap_usage = self.get_variables_of_type(RNAPUsage)
        sum_RMs = symbol_sum([x.unscaled for x in all_rnap_usage])
        free_rnap = [self.get_variables_of_type(FreeEnzyme).get_by_id(x)
                     for x in self.rnap]
        # The total RNAP capacity constraint looks like
        # ΣRMi + Σ(free RNAPj) = Σ(Total RNAPj)
        usage = sum_RMs \
                + sum([x.unscaled for x in free_rnap]) \
                - sum([x.concentration for x in self.rnap.values()])
        usage /= min([x.scaling_factor for x in self.rnap.values()])
        return usage

    def apply_rnap_catalytic_constraint(self, reaction, queue):
        """
        Given a translation reaction, apply the constraint that links it with
        RNAP usage
        :param reaction: a TranscriptionReaction
        :type reaction: TranscriptionReaction
        :return:
        """

        # Check that we indeed have a transcription reaction
        assert(isinstance(reaction, TranscriptionReaction))

        the_scaling_factor = min([x.scaling_factor for x in self.rnap.values()])

        RMi = self.add_variable(RNAPUsage,
                                reaction.gene,
                                scaling_factor = the_scaling_factor,
                                lb =0,
                                ub=1,
                                queue=False)

        # fwd_variable = reaction.forward_variable
        # bwd_variable = reaction.reverse_variable

        # v_fwd - v_bwd <= ktrans/length_aa [RNAPi]
        # v_fwd - v_bwd -  ktrans/length_aa [RNAPi] <= 0

        # Scaled in protein concentrations (eg mmol.gDW^-1)
        # σ_m is mRNA scaling factor, σ_p is protein scaling factor
        # v <= k/L [RNAP]
        # σ_m * V <= k/L * σ_m/σ_p * σ_p[RNAP]
        # v_hat <= k/L * σ_m/σ_p * [RNAP]_hat

        # v_max = self.rnap.ktrans \
        #         / reaction.nucleotide_length \
        #         * (self._mrna_scaling/self._prot_scaling) \
        #         * RMi
        # rnap_constraint_expr = fwd_variable - bwd_variable - v_max
        k = sum([this_rnap.ktrans / reaction.nucleotide_length
             for this_rnap in self.rnap.values()])

        # Scaled version:
        # rnap_constraint_expr = reaction.scaled_net - RMi # scaled
        rnap_constraint_expr = reaction.net - k * RMi.unscaled


        self.add_constraint(kind=SynthesisConstraint, hook=reaction,
                            expr=rnap_constraint_expr, ub=0,queue=queue)



    def _add_free_enzyme_ratio(self, enzyme, free_ratio):
        """
        Adds free enzyme variables to the models
        /!\ A total capacity constraint still needs to be added
        # TODO: Make that more user friendly
        :return:
        """

        # Safety_check:
        if isinstance(enzyme, Enzyme):
            enzyme = self.enzymes.get_by_id(enzyme.id)
        elif isinstance(enzyme, str):
            enzyme = self.enzymes.get_by_id(enzyme)
        else:
            raise TypeError("`enzyme' should be an enzyme object or an ID")

        # Free enzyme
        free_enz = \
            self.add_variable(kind=FreeEnzyme,
                              hook=enzyme,
                              lb=0,
                              ub=1,
                              scaling_factor = enzyme.scaling_factor)

        # Add constraint on availability of free enzymes
        expr = free_enz - free_ratio * enzyme.variable #scaled
        self.add_constraint(kind=EnzymeRatio,
                            hook=enzyme,
                            expr=expr,
                            lb=0,
                            ub=0)


    def add_ribosome(self, ribosome, free_ratio):
        """
        Adds the ribosome used by the model.

        :param ribosome:
        :type ribosome: Ribosome
        :return:
        """

        if ribosome.id in self.ribosome:
            raise KeyError('A ribosome with this id ({}) already exists in '
                           'the model'.format(ribosome.id))
        else:
            self.ribosome[ribosome.id] = ribosome

        self.add_enzymes(ribosome)

        if free_ratio > 0:
            self._add_free_enzyme_ratio(ribosome, free_ratio)

    def add_rrnas_to_rib_assembly(self, ribosome):
        """
        Adds the ribosomal RMAs to the composition of the ribosome.
        This has to be done after the transcription reactions have been added,
        so that the rRNAs synthesis reactions exist for the mass balance

        :return:
        """
        # rRNA
        rrnas = []

        f = ribosome.kdeg * ribosome.scaling_factor
        # f = self.mu_max * self.ribosome.scaling_factor

        for the_rrna_id in ribosome.rrna_composition:
            the_rrna = rRNA(id = 'rrna_' + the_rrna_id,
                            name = 'rRNA {}'.format(the_rrna_id))

            the_rrna._model = self

            rrnas.append(the_rrna)

            synthesis = self.get_transcription(the_rrna_id)

            # In the scaled mass balances, the coefficients go:
            # drRNA/dt =                                  v_tcr     - v_asm = 0
            # 1/(kdeg_rib*Erib_max) * kcat_RNAP*P_max/L * v_tcr_hat - vasm_hat = 0
            # The first coefficient is apparent stoichiometry, the second is
            # scaled automatically in the transcription reaction when rescale=True
            synthesis.add_metabolites({the_rrna:1}, rescale = True)
            # synthesis.add_metabolites({the_rrna:1}, rescale = False)


        ribosome.complexation.add_metabolites(
            {x:-1 for x in rrnas}, rescale=True)
            # {x:-1*f for x in rrnas}, rescale=False)
        self.rrnas += rrnas


    @property
    def Rt(self):
        return self.ribosome

    def _populate_ribosomes(self):
        """
        Once ribosomes have been assigned to the model, we still need to link
        them to the rest of the variables and constraints. This function creates
        the mass balance constraint on the ribosomes, as well as the total
        ribosome capacity constraint
        :return:
        """
        for this_rib in self.ribosome.values():
            self.add_rrnas_to_rib_assembly(this_rib)

            # v_complexation =   complexation.forward_variable  \
            #                  - complexation.reverse_variable

            # 1 -> Write the ribosome mass balance
            # Total amount of ribosome is in:
            # mass_balance_expr =   v_complexation            \
            #                     - self.ribosome.kdeg  * Rt  \
            #                     - self.mu             * Rt

            # Create the mass balance constraint
            # UPDATE: We now do this in add_rnap :)
            # self.add_mass_balance_constraint(this_rib.complexation,
            #                                  this_rib)

        # 2 -> Parametrize all the translation reactions with ribosomal vmax
        for trans_rxn in self.translation_reactions:
            self.apply_ribosomal_catalytic_constraint(trans_rxn)

        # 3 -> Add ribosomal capacity constraint
        self.regenerate_variables()

        ribo_usage = self._get_rib_total_capacity()
        # For nondimensionalization
        # ribo_usage = (sum_RPs + self.Rf.unscaled - self.ribosome.concentration) \
        #              / self.ribosome.scaling_factor

        # Create the capacity constraint
        self.add_constraint(kind=TotalCapacity,
                            hook=self,
                            id_='rib',
                            expr=ribo_usage,
                            lb = 0,
                            ub = 0,
                            )

        # update variable and constraints attributes
        self.regenerate_constraints()
        self.regenerate_variables()

    def _get_rib_total_capacity(self):
        free_ribosome = [self.get_variables_of_type(FreeEnzyme).get_by_id(x)
                         for x in self.ribosome]
        # CATCH : This is summing ~1500+ variable objects, and for a reason
        # sympy does not like it. Let's cut it in smaller chunks and sum
        # afterwards
        # sum_RPs = sum(self.get_variables_of_type(RibosomeUsage))
        all_ribosome_usage = self.get_variables_of_type(RibosomeUsage)
        # sum_RPs = chunk_sum(all_ribosome_usage)
        sum_RPs = symbol_sum([x.unscaled for x in all_ribosome_usage])
        # The total RNAP capacity constraint looks like
        # ΣRMi + Σ(free RNAPj) = Σ(Total RNAPj)
        ribo_usage = sum_RPs \
                     + sum([x.unscaled for x in free_ribosome]) \
                     - sum([x.concentration for x in self.ribosome.values()])
        ribo_usage /= min([x.scaling_factor for x in self.ribosome.values()])
        return ribo_usage

    def apply_ribosomal_catalytic_constraint(self, reaction):
        """
        Given a translation reaction, apply the constraint that links it with
        ribosome usage
        :param reaction: a TranslationReaction
        :type reaction: TranslationReaction
        :return:
        """

        # Check that we indeed have a translation reaction
        assert(isinstance(reaction, TranslationReaction))

        the_scaling_factor = min([x.scaling_factor for x in self.ribosome.values()])

        RPi = self.add_variable(RibosomeUsage,
                                reaction.gene,
                                scaling_factor=the_scaling_factor,
                                ub=1,
                                lb=0)

        # v_fwd - v_bwd <= kribo/length_aa [Ri]
        # v_fwd - v_bwd -  kribo/length_aa [Ri] <= 0

        # No scaling : Flux is in protein scale, and so is the ribosome concentration
        # v_max = self.ribosome.kribo \
        #         / reaction.aminoacid_length \
        #         * RPi
        # ribo_constraint_expr = fwd_variable - bwd_variable - v_max

        k = sum([rib.kribo / reaction.aminoacid_length
                 for rib in self.ribosome.values()])

        # Scaled version
        # ribo_constraint_expr = reaction.scaled_net - RPi
        ribo_constraint_expr = reaction.net - k * RPi.unscaled


        self.add_constraint(kind=SynthesisConstraint, hook=reaction,
                            expr=ribo_constraint_expr, ub=0)


    def add_genes(self, genes):
        """
        Oddly I could not find this method in cobra. Adds one or several genes
        to the model.

        :param genes:
        :type genes: Iterable(Gene) or Gene
        :return:
        """
        if not hasattr(genes,'__iter__'):
            genes = [genes]

        new_genes = [x for x in genes if x.id not in self.genes]
        modified_genes = [x for x in genes if x.id in self.genes]

        for g in new_genes:
            self._add_gene(g)

        for g in modified_genes:
            raise NotImplementedError()

    def _add_gene(self, gene):
        gene._model = self
        self.genes.append(gene)


    #-------------------------------------------------------------------------#

    def sanitize_varnames(self):
        """
        Makes variable name safe for the solvers. In particular, variables whose
        name start with
        :return:
        """
        for met in self.metabolites:
            if met.id[0].isdigit():
                met.id = '_'+met.id
                self.logger.info('Sanitized variable name {}'.format(met.id))
        for rxn in self.reactions:
            if rxn.id[0].isdigit():
                rxn.id = '_'+rxn.id
                self.logger.info('Sanitized variable name {}'.format(rxn.id))

        Model.repair(self)

    def print_info(self, specific = False):
        """
        Print information and counts for the cobra_model
        :return:
        """
        if not specific:
            LCSBModel.print_info(self)

        n_reactions = len(self.reactions)
        n_enzymes = len(self.enzymes)
        n_enzymatic_reactions   = len([x for x in self.reactions \
                                       if isinstance(x, EnzymaticReaction)])

        info = pd.DataFrame(columns = ['value'])
        info.loc['num enzymes'] = n_enzymes
        info.loc['num enzymatic_reactions'] = n_enzymatic_reactions
        info.loc['pct enzymatic_reactions'] = n_enzymatic_reactions/n_reactions*100
        info.index.name = 'key'

        print(info)

    def __deepcopy__(self, memo):
        """
        Calls self.copy() to return an independant copy of the model

        :param memo:
        :return:
        """

        return self.copy()

    def copy(self):
        """
        Pseudo-smart copy of the model using dict serialization. This builds a
        new model from the ground up, with independwnt variables, solver, etc.

        :return:
        """

        from ..io.dict import model_from_dict, model_to_dict
        dictmodel = model_to_dict(self)
        new = model_from_dict(dictmodel)

        copy_solver_configuration(self, new)

        return new
