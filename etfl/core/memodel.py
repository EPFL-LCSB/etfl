# -*- coding:utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: ETFL team

Core for the ME-part

"""

import numpy as np
import optlang
import pandas as pd
import sympy
from cobra import Model, Reaction, Gene
from cobra.core import Solution, DictList
from collections import defaultdict
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
from .expression import build_trna_charging, \
    make_stoich_from_aa_sequence, make_stoich_from_nt_sequence, \
    degrade_peptide, degrade_mrna
from ..optim.constraints import CatalyticConstraint, ForwardCatalyticConstraint,\
    BackwardCatalyticConstraint, EnzymeMassBalance, \
    rRNAMassBalance, mRNAMassBalance, tRNAMassBalance, DNAMassBalance, \
    GrowthCoupling, TotalCapacity, ExpressionCoupling, EnzymeRatio, \
    GrowthChoice, EnzymeDegradation, mRNADegradation,\
    LinearizationConstraint, SynthesisConstraint, SOS1Constraint,\
    InterpolationConstraint
from ..optim.variables import ModelVariable, GrowthActivation, \
    EnzymeVariable, LinearizationVariable, RibosomeUsage, RNAPUsage, \
    FreeRibosomes, BinaryActivator, InterpolationVariable, DNAVariable, \
    GrowthRate, GenericVariable

from pytfa.core.model import LCSBModel
from pytfa.optim.reformulation import petersen_linearization
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
        self.enzymes = DictList()
        self.mrnas = DictList()
        self.rrnas = DictList()
        self.trnas = DictList()
        self.peptides = DictList()
        self.transcription_reactions = DictList()
        self.translation_reactions = DictList()
        self.complexation_reactions = DictList()
        self.degradation_reactions = DictList()


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
        mu in [0.1, 0.9] with nbins = 8
        => mu = 0.15 OR mu = 0.25 OR ... OR mu = 0.85

        Using binary exapnsion of the bins instead of a list of 0-1s
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

            free_pep = Peptide( id=gene_id,
                                name='Peptide, {}'.format(gene_id),
                                gene_id=gene_id)

            free_pep._model = self
            self.peptides += [free_pep]



    def add_dummies(self, nt_ratios, mrna_kdeg, mrna_length, aa_ratios,
                    enzyme_kdeg, peptide_length):
        """

        create dummies to enforce mrna and peptide production even in the
        absence of data for all mrnas and proteins
        absence of data for all mrnas and proteins


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

        gtp = self.essentials['gtp']
        gdp = self.essentials['gdp']
        h2o = self.essentials['h2o']
        h = self.essentials['h']
        ppi = self.essentials['ppi']

        # Create a dummy gene and override the sequences with input data
        dummy_gene = ExpressedGene(id='dummy_gene',
                                   name='Dummy Gene',
                                   sequence='')
        dummy_gene._rna = 'N'*mrna_length
        dummy_gene._peptide = 'X'*peptide_length

        self.add_genes([dummy_gene])

        # Create a dummy mRNA
        dummy_mrna = mRNA(id='dummy_gene',
                          name='dummy mRNA',
                          kdeg=mrna_kdeg,
                          gene_id = dummy_gene.id)

        nt_weights = [v*molecular_weight(k, 'RNA') for k,v in nt_ratios.items()]
        dummy_mrna.molecular_weight = mrna_length*sum(nt_weights)/ 1000 # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1

        self.add_mrnas([dummy_mrna], add_degradation = False)

        dummy_transcription = TranscriptionReaction(id=self._get_transcription_name(dummy_mrna.id),
                                                    name = 'Dummy Transcription',
                                                    gene_id=dummy_gene.id,
                                                    enzymes=self.rnap,
                                                    scaled=True)
        self.add_reactions([dummy_transcription])
        self.transcription_reactions += [dummy_transcription]

        # Use the input ratios to make the stoichiometry
        transcription_mets = {
            self.metabolites.get_by_id(self.rna_nucleotides[k]):
                -1 * v * mrna_length
            for k,v in nt_ratios.items()
        }
        transcription_mets[ppi] = mrna_length

        dummy_transcription.add_metabolites(transcription_mets, rescale = True)


        # Add the degradation
        mrna_deg_stoich = {
            self.metabolites.get_by_id(self.rna_nucleotides_mp[k]):
                -1 * v * mrna_length
            for k,v in nt_ratios.items()
        }
        mrna_deg_stoich[h2o] = -1 * mrna_length
        mrna_deg_stoich[h] = 1 * mrna_length

        self._make_degradation_reaction(deg_stoich=mrna_deg_stoich,
                                        macromolecule=dummy_mrna,
                                        kind=mRNADegradation,
                                        scaled = True)

        self.add_mass_balance_constraint(dummy_transcription, dummy_mrna)

        # Create a dummy peptide
        dummy_peptide = Peptide(id='dummy_peptide',
                                name='Dummy peptide',
                                gene_id=dummy_gene.id)

        aa_weights = [v*molecular_weight(k, 'protein') for k,v in aa_ratios.items()]

        dummy_peptide.molecular_weight = peptide_length*sum(aa_weights)/ 1000 # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1


        # Create a dummy protein made of the dummy peptide
        dummy_protein = Enzyme(id='dummy_enzyme',
                               name='Dummy Enzyme',
                               kcat=0,
                               composition=[dummy_peptide.id],
                               kdeg=enzyme_kdeg)

        self.add_enzymes([dummy_protein])
        dummy_peptide._model = self
        self.peptides += [dummy_peptide]

        # Making the reactions and adding them to the model so that they can
        # reference genes

        dummy_translation = TranslationReaction(id=self._get_translation_name(dummy_peptide),
                                                name='Dummy Translation',
                                                gene_id=dummy_gene.id,
                                                enzymes=self.ribosome,
                                                scaled=True)
        dummy_complexation = ProteinComplexation(id='dummy_complexation',
                                                 name='Dummy Complexation',
                                                 target=dummy_protein,
                                                 scaled=True)
        self.add_reactions([dummy_translation, dummy_complexation])
        self.translation_reactions += [dummy_translation]
        self.complexation_reactions += [dummy_complexation]

        # Use the input ratios to make the stoichiometry
        translation_mets = {}

        for k,v in aa_ratios.items():
            the_met_id = self.aa_dict[k]
            the_charged_trna, the_uncharged_trna, _ = self.trna_dict[the_met_id]
            translation_mets[the_charged_trna  ] = -1*v*peptide_length
            translation_mets[the_uncharged_trna] =  1*v*peptide_length


        translation_mets[self.metabolites.get_by_id(gtp)] = -2*peptide_length
        translation_mets[self.metabolites.get_by_id(h2o)] = -2*peptide_length
        translation_mets[self.metabolites.get_by_id(gdp)] =  2*peptide_length
        translation_mets[self.metabolites.get_by_id( h )] =  2*peptide_length

        # Do not forget to extract the tRNAs from the stoichiometry, since they
        # get diluted
        self._extract_trna_from_reaction(translation_mets, dummy_translation)

        dummy_translation.add_metabolites(translation_mets, rescale=True)
        dummy_translation.add_peptide(dummy_peptide)

        dummy_complexation.add_peptides({dummy_peptide:-1})
        dummy_protein.init_variable()

        # Finally add the degradation flux

        prot_deg_stoich = dict()

        for k, v in aa_ratios.items():
            the_met_id = self.aa_dict[k]
            prot_deg_stoich[the_met_id] = v * peptide_length
        prot_deg_stoich[h2o] = -1*peptide_length

        self._make_degradation_reaction(deg_stoich=prot_deg_stoich,
                                        macromolecule=dummy_protein,
                                        kind=EnzymeDegradation,
                                        scaled=True)

        # Now we can add the mass balance constraint

        self.add_mass_balance_constraint(dummy_complexation, dummy_protein)

    def add_interpolation_variables(self):
        lambdas = []
        for e in range(self.n_mu_bins):
            lambda_i = self.add_variable(kind=BinaryActivator,
                                         hook=self,
                                         id_=str(e),
                                         lb=0,
                                         ub=1
                                         )
            lambdas += [lambda_i]
        sos_expr = symbol_sum(lambdas)

        self.add_constraint(kind=SOS1Constraint,
                            hook=self,
                            id_='interpolation_integer_SOS1',
                            expr=sos_expr,
                            lb=1,
                            ub=1)

        ga_vars = self.get_ordered_ga_vars()
        # mu_integer is the fraction coefficient of mu/mu_max:
        # mu_integer = delta_0*2^0 + delta_1*2^1 + ... + delta_n*2^n
        the_mu_integer = symbol_sum([(2 ** i) * ga_i
                                     for i, ga_i in enumerate(ga_vars)])

        # We want to equate the mu_integer with the currently active lambda index
        # 0*lambda_0 + 1*lambda_1 + ... + n*lambda_n = mu_integer
        ic_expr = symbol_sum([e*l for e,l in enumerate(lambdas)]) - the_mu_integer

        self.add_constraint(kind=InterpolationConstraint,
                            hook=self,
                            id_='growth_activators__EQ__interpolation_integers',
                            expr=ic_expr,
                            lb=0,
                            ub=0)

        self.regenerate_variables()
        self.regenerate_constraints()


    def add_protein_mass_requirement(self, mu_values, p_rel):
        """
        Adds protein synthesis requirement

        input of type:
        mu_values=[ 0.6,        1.0,        1.5,        2.0,        2.5     ]
        p_rel   = [ 0.675676,   0.604651,   0.540416,   0.530421,   0.520231]

        mu_values in [h^-]
        p_rel in [g/gDw]

        :param mu_values:
        :param p_rel:
        :return:
        """

        activation_vars = self.get_variables_of_type(BinaryActivator)

        model_mus = [x[0] for x in self.mu_bins]
        p_hat = np.interp(x= model_mus,
                          xp=mu_values,
                          fp=p_rel)

        p_ref = symbol_sum([x*y for x,y in zip(p_hat, activation_vars)])

        # mmol.gDw^-1 / [scaling]
        enzyme_vars    = self.enzymes.list_attr('concentration')
        # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
        enzyme_weights = self.enzymes.list_attr('molecular_weight')

        tot_prot = symbol_sum([x*y for x,y in zip(enzyme_weights,enzyme_vars)])

        # For legibility
        prot_ggdw = self.add_variable(kind=InterpolationVariable, hook=self,
                                      id_='prot_ggdw',
                                      lb=0,
                                      ub=1, # can't have more prot than cell mass
                                      )

        # MW_1*[E1] + MW_2*[E2] + ... + MW_n*[En] = prot_ggdw
        mass_variable_def = tot_prot - prot_ggdw

        # E_ggdw = E_ref
        mass_coupling_expr = prot_ggdw - p_ref

        epsilon = max(abs(np.diff(p_hat)))

        self.add_constraint(kind=InterpolationConstraint,
                            hook=self,
                            id_='prot_weight_definition',
                            expr=mass_variable_def,
                            lb=0,
                            ub=0,
                            )

        self.add_constraint(kind=InterpolationConstraint,
                            hook=self,
                            id_='prot_interpolation',
                            expr=mass_coupling_expr,
                            lb=-1 * epsilon,
                            ub=epsilon,
                            )

        self.interpolation_protein = p_hat
        self._interpolation_protein_tolerance = epsilon

        self.regenerate_variables()
        self.regenerate_constraints()


    def add_rna_mass_requirement(self, mu_values, rna_rel):
        """
        Adds RNA synthesis requirement

        input of type:
        mu_values = [   0.6,        1.0,        1.5,        2.0,        2.5     ]
        rna_rel   = [   0.135135    0.151163    0.177829    0.205928    0.243931]

        mu_values in [h^-]
        rna_rel in [g/gDw]

        :param mu_values:
        :param rna_rel:
        :return:
        """

        activation_vars = self.get_variables_of_type(BinaryActivator)

        model_mus = [x[0] for x in self.mu_bins]
        m_hat = np.interp(x= model_mus,
                          xp=mu_values,
                          fp=rna_rel)

        m_ref = symbol_sum([x*y for x,y in zip(m_hat, activation_vars)])

        rna_vars    = self.mrnas.list_attr('concentration') # mmol.gDw^-1 / [scaling]
        rna_weights = self.mrnas.list_attr('molecular_weight') # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1

        tot_rna = symbol_sum([x*y for x,y in zip(rna_weights,rna_vars)])


        # For legibility
        mrna_ggdw = self.add_variable(kind=InterpolationVariable,
                                      hook=self,
                                      id_='mrna_ggdw',
                                      lb=0,
                                      ub=1, #can't have more rna than cell mass
                                      )

        # MW_1*[rna1] + MW_2*[rna2] + ... + MW_n*[rna_n] = mRNA_ggdw
        mass_variable_def = tot_rna - mrna_ggdw

        # mRNA_ggdw = mRNA_ref
        mass_coupling_expr = mrna_ggdw - m_ref

        epsilon = max(abs(np.diff(m_hat)))

        self.add_constraint(kind=InterpolationConstraint,
                            hook=self,
                            id_='mRNA_weight_definition',
                            expr=mass_variable_def,
                            lb=0,
                            ub=0,
                            )

        self.add_constraint(kind=InterpolationConstraint,
                            hook=self,
                            id_='mRNA_interpolation',
                            expr=mass_coupling_expr,
                            lb=-1*epsilon,
                            ub=epsilon,
                            )


        self.interpolation_mrna = m_hat
        self._interpolation_mrna_tolerance = epsilon

        self.regenerate_variables()
        self.regenerate_constraints()


    def add_dna_mass_requirement(self, mu_values, dna_rel, gc_ratio,
                                 chromosome_len, dna_dict, ppi='ppi_c'):
        """
        Adds DNA synthesis requirement

        input of type:
        mu_values = [   0.6,        1.0,        1.5,        2.0,        2.5     ]
        dna_rel   = [   0.135135    0.151163    0.177829    0.205928    0.243931]

        mu_values in [h^-]
        dna_rel in [g/gDw]

        :param mu_values:
        :param dna_rel:
        :return:
        """

        self.dna_nucleotides = dna_dict

        # Get mu interpolation
        activation_vars = self.get_variables_of_type(BinaryActivator)

        model_mus = [x[0] for x in self.mu_bins]
        m_hat = np.interp(x= model_mus,
                          xp=mu_values,
                          fp=dna_rel)

        m_ref = symbol_sum([x*y for x,y in zip(m_hat, activation_vars)])


        # Add DNA variable:

        dna = DNA(kdeg = 0, dna_len=chromosome_len, gc_ratio=gc_ratio)
        # Assumption: kdeg for DNA is close to 0
        self.add_dna(dna)

        # Create dummy DNA reaction
        dna_formation = DNAFormation(id='DNA_formation', name = 'DNA Formation',
                                     dna=dna, scaled = True)
        self.add_reactions([dna_formation])

        # In this formulation, we make 1 unit of the whole chromosome with NTPs
        g = gc_ratio
        mets = {v: -1 * chromosome_len * (g if k.lower() in 'gc' else 1 - g)
                for k, v in self.dna_nucleotides.items()}
        # Don't forget to release ppi (2 ppi per bp)
        mets[ppi] = 2 * chromosome_len


        dna_formation.add_metabolites(mets)

        # Add mass balance : 0 = v_syn - [mu]*[DNA]
        self.add_mass_balance_constraint(
            synthesis_flux=dna_formation,
            macromolecule=dna)

        # DNA mass (BioPython has g.mol^-1, while we are in mmol)
        ma = molecular_weight('A', seq_type='DNA')/ 1000 # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
        mt = molecular_weight('T', seq_type='DNA')/ 1000 # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
        mc = molecular_weight('C', seq_type='DNA')/ 1000 # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
        mg = molecular_weight('G', seq_type='DNA')/ 1000 # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1

        #              g.mmol(bp)^-1        * mmol(bp)/mmol(dna) * mmol(dna).gDW^-1
        tot_dna = ( (1-g)*(ma+mt) + g*(mc+mg) ) * chromosome_len * dna.concentration

        # For legibility
        dna_ggdw = self.add_variable(kind=InterpolationVariable,
                                     hook=self,
                                     id_='dna_ggdw',
                                     lb=0,
                                     ub=1, #can't have more dna than cell mass
                                     )

        # MW_avg*[DNA] = mRNA_ggdw
        # 1/scaling because the [X]s are scaled (eg mmol.ggDW^-1 -> back to mol.ggDW^1)
        mass_variable_def  = tot_dna - dna_ggdw

        # DNA_ggdw = DNA_ref
        mass_coupling_expr = dna_ggdw - m_ref

        epsilon = max(abs(np.diff(m_hat)))

        self.add_constraint(kind=InterpolationConstraint,
                            hook=self,
                            id_='DNA_weight_definition',
                            expr=mass_variable_def,
                            lb=0,
                            ub=0,
                            )

        self.add_constraint(kind=InterpolationConstraint,
                            hook=self,
                            id_='DNA_interpolation',
                            expr=mass_coupling_expr,
                            lb=-1*epsilon,
                            ub=epsilon,
                            )


        self.interpolation_dna = m_hat
        self._interpolation_dna_tolerance = epsilon

        self.regenerate_variables()
        self.regenerate_constraints()

    def add_essentials(self, essentials, aa_dict, rna_nucleotides,
                       rna_nucleotides_mp):
        """
        Marks important metabolites for expression

        :param essentials:

        :param aa_dict: A dictionnary of aminoacid letter to amicoacid met id
            Example :
            ```python
            aa_dict = {
                        'A':'ala__L_c',
                        'R':'arg__L_c',
                        ...
                    }
            ```
        :param rna_nucleotides: A dictionnary of RNA nucleotide letter to nucleotideTP met id
            Example :
            ```python
            rna_nucleotides = {
                        'A':'atp_c',
                        'U':'utp_c',
                        ...
                    }
            ```
        :param rna_nucleotides_mp:
        :return:
        """

        self.essentials = essentials
        self.aa_dict = aa_dict
        self.rna_nucleotides = rna_nucleotides
        self.rna_nucleotides_mp = rna_nucleotides_mp


    def build_expression(self):
        """
        Given a dictionnary from amino acids nucleotides to metabolite names,
        goes through the list of genes in the model that have sequence
        information to build transcription and traduction reactions
        :return:
        """

        aa_dict = self.aa_dict
        rna_nucleotides = self.rna_nucleotides

        atp = self.essentials['atp']
        amp = self.essentials['amp']
        ppi = self.essentials['gtp']
        h2o = self.essentials['h2o']
        h = self.essentials['h']

        self.trna_dict = build_trna_charging(self,aa_dict,atp,amp,ppi,h2o,h)
        self.add_trnas([item for sublist in self.trna_dict.values()
                        for item in sublist if isinstance(item, tRNA)])

        # Check that the ribosomes have been added
        if self.ribosome is None:
            raise Exception(
                'A ribosome has to be added with the add_ribosome method')

        # Check that the RNAP has been added
        if self.rnap is None:
            raise Exception(
                'A RNA Polymerase has to be added with the add_rnap method')

        self.express_genes(self.genes)

    def express_genes(self, gene_list):

        for gene in gene_list:

            if isinstance(gene, str):
                the_gene = self.genes.get_by_id(gene)
            else:
                the_gene = gene

            if not isinstance(the_gene, ExpressedGene):
                continue

            # Build the translation
            self._add_gene_translation_reaction(the_gene)
            # Build the transcription
            self._add_gene_transcription_reaction(the_gene)

    def _add_gene_translation_reaction(self, gene):
        """

        :param gene: A gene of the model that has sequence data
        :type gene: etfl.core.ExpressedGene
        :return:
        """

        gtp = self.essentials['gtp']
        gdp = self.essentials['gdp']
        pi  = self.essentials['pi']
        h2o = self.essentials['h2o']
        h   = self.essentials['h']

        Rmax = 1 / self.ribosome.molecular_weight

        rxn = TranslationReaction(
            id='{}_translation'.format(gene.id),
            name='Translation, {}'.format(gene.id),
            gene_id= gene.id,
            enzymes=self.ribosome,
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

        self._extract_trna_from_reaction(aa_stoichiometry, rxn)

        rxn.add_metabolites(aa_stoichiometry, rescale=True)

        free_peptide = self.peptides.get_by_id(gene.id)

        rxn.add_peptide(free_peptide)

        # Add ribosome as necessary enzyme
        rxn.gene_reaction_rule = self.ribosome.id
        self.translation_reactions += [rxn]

    def _extract_trna_from_reaction(self, aa_stoichiometry, rxn):
        # Extract the tRNAs, since they will be used for a different mass balance
        # in self.add_trna_mass_balances
        for met, stoich in list(aa_stoichiometry.items()):
            if isinstance(met, tRNA):
                rxn.trna_stoich[met.id] = aa_stoichiometry.pop(met)

    def _add_gene_transcription_reaction(self, gene):
        """

        :param gene: A gene of the model that has sequence data
        :type gene: etfl.core.ExpressedGene
        :return:
        """

        ppi  = self.essentials['ppi']

        nt_stoichiometry = make_stoich_from_nt_sequence(gene.rna,
                                                        self.rna_nucleotides,
                                                        ppi
                                                        )

        rxn = TranscriptionReaction(
            id=self._get_transcription_name(gene.id),
            name='Transcription, {}'.format(gene.id),
            gene_id= gene.id,
            enzymes=self.rnap,
            upper_bound=1,
            scaled=True)
        self.add_reactions([rxn])

        rxn.add_metabolites(nt_stoichiometry)

        # Add rnap as necessary enzyme
        rxn.gene_reaction_rule = ' & '.join(self.rnap.composition)


        self.transcription_reactions += [rxn]

    def add_trna_mass_balances(self):
        """
        Once the tRNAs, transcription and translation reactions have been added,
        we need to add the constraints:

        d/dt [charged_tRNA]   =  v_charging - sum(nu_trans*v_trans) - mu*[charged_tRNA]
        d/dt [uncharged_tRNA] = -v_charging + sum(nu_trans*v_trans) - mu*[uncharged_tRNA]

        The stoichiometries are set from the reaction dict in _extract_trna_from_reaction

        We also need to scale the tRNAs in mRNA space and unscale the translation:

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
        self.coupling_dict = coupling_dict
        self.add_enzymes(coupling_dict.values())

        # /!\ We modify the reaction list
        # self.add_gene_reactions()

        # Generic reactions <-> Enzymes coupling
        for rid in tqdm(self.coupling_dict, desc='cat. constraints'):
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
                self.logger.error('Could not find reaction {} in the coupling dictionnary'.format(r.id))

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

        complexation = self.add_complexation_from_enzymes(reaction.enzymes)

        # Why is that happening ?
        if not complexation:
            return None

        v_max_fwd = dict()
        v_max_bwd = dict()

        # Write v_max constraint
        fwd_variable = reaction.forward_variable
        bwd_variable = reaction.reverse_variable

        zipped_prot_isoz = self.match_enzymes_to_complexes(reaction.enzymes,
                                                               complexation)

        protein2isozyme_dict = dict(zipped_prot_isoz)

        for e, (enz, comp) in enumerate(protein2isozyme_dict.items()):
            # If the enzymes has the same kcat for both directions
            # v_fwd  <= kcat_fwd [E]
            # v_fwd - kcat_fwd [E] <= 0

            v_max_fwd[e] = enz.kcat_fwd * enz.concentration
            v_max_bwd[e] = enz.kcat_bwd * enz.concentration

            self.add_mass_balance_constraint(comp, enz, queue=False)

            comp.enzyme = enz
            enz.complexation = comp

        k_f = max([x.kcat_fwd for x in protein2isozyme_dict])
        k_b = max([x.kcat_bwd for x in protein2isozyme_dict])
        # k_f = np.median([x.kcat_fwd for x in self.enzymes])
        # k_b = np.median([x.kcat_bwd for x in self.enzymes])
        E_m = max([x.scaling_factor for x in protein2isozyme_dict])

        # v_fwd <= sum(kcat_i*E_i)
        # for all i, E_i <= E_max (= 1g/gDW)
        # v_fwd / sum(kcat_i*E_i^max) <= sum(kcat_i*E_i) / sum(kcat_i*E_i^max) (<= 1)

        enz_constraint_expr_fwd = (fwd_variable - sum(v_max_fwd.values()))/(k_f*E_m)
        enz_constraint_expr_bwd = (bwd_variable - sum(v_max_bwd.values()))/(k_b*E_m)
        enz_constraint_expr_bwd = (bwd_variable - sum(v_max_bwd.values()))/(k_b*E_m)


        self.add_constraint(kind=ForwardCatalyticConstraint, hook=reaction,
                            expr=enz_constraint_expr_fwd, ub=0, queue=True)
        self.add_constraint(kind=BackwardCatalyticConstraint, hook=reaction,
                            expr=enz_constraint_expr_bwd, ub=0, queue=True)

    def add_mass_balance_constraint(self, synthesis_flux, macromolecule=None,
                                    queue=False):
        """
        Adds a mass balance constraint of the type
        d[E]/dt = 0 <=> v_synthesis - k_deg*[M] - μ*[M] = 0
        for a macromolecule (mRNA or enzyme)
        :param synthesis_flux:
        :param macromolecule:
        :return:
        """

        kwargs = dict()

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
            mu_ub = self.mu_max
            sigma = mu_ub * macromolecule.scaling_factor
            # The synthesis flux comes unscaled [flux units]
            # expr = 1/sigma * synthesis_flux - 1/mu_ub * z
            expr = synthesis_flux - macromolecule.scaling_factor * z
        elif isinstance(macromolecule, DNA):
            kind = DNAMassBalance
            kwargs['id_'] = 'dna'
            hook = self
            mu_ub = self.mu_max
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

        # z = lambda_i * E_hat <= 1
        ub = 1

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
            z_name = '__MUL__'.join([ga_i.name, E.name])
            # Add the variables
            model_z_i = self.add_variable(kind=LinearizationVariable,
                                          hook=self,
                                          id_=z_name,
                                          lb=0,
                                          ub=ub,
                                          queue=False)

            # z_i, cons = glovers_linearization(b = ga_i, fy=E, L=E.lb, U=E.ub, z=model_z_i)
            z_i, new_constraints = petersen_linearization(b=ga_i, x=E, M=E.ub,
                                                          z=model_z_i)

            # Add the constraints:
            for cons in new_constraints:
                # Do not forget to substitute the sympy symbol in the constraint
                # with a variable  !
                # new_expression = cons.expression.subs(z_i, model_z_i.variable)
                # EDIT: Not anymore needed if we supply the variable

                self.add_constraint(kind=LinearizationConstraint,
                                    hook=self,
                                    id_=cons.name,
                                    expr=cons.expression,
                                    # expr=new_expression,
                                    ub=cons.ub,
                                    lb=cons.lb,
                                    queue=queue)

            out_expr += (2 ** i) * self.mu_approx_resolution * model_z_i

        return out_expr

    def get_ordered_ga_vars(self):
        # ga_i is a binary variable for the binary expansion f the fraction on N
        # of the max growth rate
        ga_vars = self.get_variables_of_type(GrowthActivation)
        ga_vars = sorted(ga_vars, key=lambda x: x.ix)
        return ga_vars


    def add_complexation_from_enzymes(self,enzymes):
        """
        Reads Enzyme.composition to find complexation reaction from enzyme information

        :param reaction:
        :type reaction: cobra.Reaction
        :return:
        """

        complexation = []

        for e,this_isozyme in enumerate(enzymes):
            # Does this enzyme already have a complexation reaction ?
            # This happens if an enzyme is used in several reactions
            if this_isozyme.complexation is not None:
                complexation += [this_isozyme.complexation]
            else:
                if not this_isozyme.composition:
                    self.logger.warning('Enzyme {} has no composition'
                                        .format(this_isozyme.id))
                    continue

                this_id = '{}_complex_{}'.format(this_isozyme.id,e)
                this_name = '{} Complexation {}'.format(this_isozyme.id,e)


                this_complexation = ProteinComplexation(id = this_id,
                                                        name = this_name,
                                                        target=this_isozyme,
                                                        # upper_bound=1,
                                                        scaled=True)

                peptides = {self.peptides.get_by_id(k):-v \
                            for k,v in this_isozyme.composition.items()}

                self.add_reactions([this_complexation])
                this_complexation.add_peptides(peptides)

                complexation += [this_complexation]
                this_isozyme.init_variable()

                # also add degradation
                self._add_enzyme_degradation(this_isozyme, scaled=True, queue=True)

                # Add it to a specific index
                self.complexation_reactions += [this_complexation]
        # self._push_queue()

        return complexation

    def match_enzymes_to_complexes(self, enzymes, complexes):
        #TODO: Implement this better
        # if only one prot, replicate it with similar kdeg and kcat
        if   len(enzymes) == len(complexes):
            return zip(enzymes, complexes)
        elif len(enzymes) == 1:
            enzyme_list = self.replicate_enzyme(enzymes[0], len(complexes))
            return zip(enzyme_list, complexes)
        else:
            raise NotImplementedError


    def add_enzymes(self, enzyme_list):
        """
        Adds an Enzyme object, or iterable of Enzyme objects, to the model
        :param enzyme_list:
        :type enzyme_list:Iterable(Enzyme) or Enzyme
        :return:
        """
        if not hasattr(enzyme_list, '__iter__'):
            enzyme_list = [enzyme_list]
        else:
            enzyme_list = list(enzyme_list)
        if len(enzyme_list) == 0:
            return None

        if not isinstance(enzyme_list[0],Enzyme):
            enzyme_list = [x for item in enzyme_list for x in item]

        # First check whether the enzymes exist in the model
        # Also enzymes could be declared twice for different reactions,
        # hence turn the list into a set
        enzyme_list = [x for x in set(enzyme_list) if x.id not in self.enzymes]

        for enz in enzyme_list:
            enz._model = self

        self.enzymes += enzyme_list


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
            self.enzymes.pop(enz.id)

    def replicate_enzyme(self, enzyme, n_replicates):
        """
        Replicates an enzyme n_replicates times, with similar kcat and kdeg.
        Useful for isozymes

        :param enzyme:
        :type enzyme: pytfa.me.Enzyme
        :param n_replicates:
        :type n_replicates: int
        :return:
        """
        self.remove_enzymes([enzyme])

        new_enzymes = list()
        for e in range(n_replicates):
            new_enz = Enzyme(id =enzyme.id + '_{}'.format(e),
                             kcat_fwd = enzyme.kcat_fwd,
                             kcat_bwd = enzyme.kcat_bwd,
                             kdeg = enzyme.kdeg,
                             # name = enzyme.name + ' - Replicate {}'.format(e)
                             )
            new_enzymes.append(new_enz)

        self.add_enzymes(new_enzymes)
        return new_enzymes

    def _add_enzyme_degradation(self, enzyme, scaled=True, queue=False):

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

        h2o = self.essentials['h2o']
        h = self.essentials['h']

        if mrna.kdeg is None or np.isnan(mrna.kdeg):
            return None

        degradation_mets = degrade_mrna(mrna, self.rna_nucleotides_mp, h2o, h)

        self._make_degradation_reaction(degradation_mets,mrna,mRNADegradation,
                                        scaled=scaled, queue=queue)


    def _make_degradation_reaction(self, deg_stoich, macromolecule,
                                   kind, scaled, queue=False):

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


    def populate_expression(self):
        """
        Add the coupling between mRNA availability and ribosome charging
        The number of ribosomes assigned to a mRNA species is lower than
        the number of such mRNA times the max number of ribosomes that can sit
        on the mRNA:
        [RPi] <= loadmax_i*[mRNAi]

        loadmax is : len(peptide_chain)/occupation(ribo)
        "Their distance from one another along the mRNA is at least the size
        of the physical footprint of a ribosome (≈20 nm, BNID 102320, 100121)
        which is the length of about 60 base pairs (length of
        nucleotide ≈0.3 nm, BNID 103777), equivalent to ≈20 aa."
        "http://book.bionumbers.org/how-many-proteins-are-made-per-mrna-molecule/"

        hence:
        [RPi] <= L_nt/Ribo_footprint * [mRNA]

        :return:
        """
        self._populate_rnap()
        self._populate_ribosomes()

        ribo_footprint_size = 60 # see docstring

        self._push_queue()

        for the_mrna in tqdm(self.mrnas, desc='populating expression'):

            # Get the synthesis_flux
            syn_id = self._get_transcription_name(the_mrna.id)
            syn = self.transcription_reactions.get_by_id(syn_id)
            # Add the mass balance constraint for the mrna
            self.add_mass_balance_constraint(syn, the_mrna, queue=True)

            # Get the ribosomes assigned to this translation
            RPi_hat = getattr(self, camel2underscores(RibosomeUsage.__name__)) \
                .get_by_id(the_mrna.id)

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
            # expression_coupling = RPi_hat \
            #                       - polysome_size * scaling_factor * mrna_hat
            expression_coupling = RPi_hat - polysome_size * mrna_hat


            # Add expression coupling
            self.add_constraint(kind = ExpressionCoupling,
                                hook = the_mrna,
                                expr = expression_coupling,
                                queue = True,
                                ub = 0)

        self._push_queue()
        self.regenerate_variables()
        self.regenerate_constraints()

    def _get_transcription_name(self, the_mrna_id):
        return '{}_transcription'.format(the_mrna_id)

    def _get_translation_name(self, the_peptide_id):
        return '{}_translation'.format(the_peptide_id)

    def get_translation(self, the_peptide_id):
        return self.reactions.get_by_id(self._get_translation_name(the_peptide_id))

    def get_transcription(self, the_peptide_id):
        return self.reactions.get_by_id(self._get_transcription_name(the_peptide_id))

    def add_rnap(self, rnap, free_ratio = 0.8):
        """
        Adds the RNA Polymerase used by the model.

        :param rnap:
        :type rnap: pytfa.me.Ribosome
        :return:
        """

        self.rnap = rnap

        self.add_enzymes(rnap)

        # 0 -> We still need to add the virtual complexation of RNA Polymerases:
        peptide_stoich = defaultdict(int)
        for rprot_id in rnap.composition:
            peptide_stoich[self.peptides.get_by_id(rprot_id)] -= 1

        complexation = ProteinComplexation(id='rnap_complex',
                                           name='RNA Polymerase complexation',
                                           target=rnap,
                                           scaled=True)

        self.add_reactions([complexation])
        complexation.add_peptides(peptide_stoich)

        self.complexation_reactions += [complexation]
        self.rnap.init_variable()

        # Let us not forget the degradation
        self._add_enzyme_degradation(rnap, scaled=True)
        self.init_rnap_variables(free_ratio=free_ratio)


    def _populate_rnap(self):
        """
        Once ribosomes have been assigned to the model, we still need to link
        them to the rest of the variables and constraints. This function creates
        the mass balance constraint on the ribosomes, as well as the total
        ribosome capacity constraint
        :return:
        """

        # v_complexation =   complexation.forward_variable  \
        #                  - complexation.reverse_variable

        # 1 -> Write the RNAP mass balance
        # Create the mass ba;ance constraint
        self.add_mass_balance_constraint(self.rnap.complexation, self.rnap)

        # 2 -> Parametrize all the transcription reactions with RNAP vmax
        for trans_rxn in self.transcription_reactions:
            self.apply_rnap_catalytic_constraint(trans_rxn)

        # 3 -> Add RNAP capacity constraint
        self.regenerate_variables()

        all_rnap_usage = self.get_variables_of_type(RNAPUsage)
        sum_RMs = symbol_sum([x.unscaled for x in all_rnap_usage])

        # For nondimensionalization, we divide by the sum of D upperbounds
        D = self.rnap.scaling_factor
        # usage = (sum_RMs- self.rnap.concentration) / D
        usage = (sum_RMs + self.RNAPf.unscaled - self.rnap.concentration) / \
                self.rnap.scaling_factor

        # Create the capacity constraint
        self.add_constraint(kind=TotalCapacity,
                            hook=self.rnap,
                            expr=usage,
                            lb = 0,
                            ub = 0,
                            )

        # update variable and constraints attributes
        self.regenerate_constraints()
        self.regenerate_variables()

    def apply_rnap_catalytic_constraint(self, reaction):
        """
        Given a translation reaction, apply the constraint that links it with
        RNAP usage
        :param reaction: a TranscriptionReaction
        :type reaction: TranscriptionReaction
        :return:
        """

        # Check that we indeed have a transcription reaction
        assert(isinstance(reaction, TranscriptionReaction))

        RMi = self.add_variable(RNAPUsage,
                                reaction.gene,
                                scaling_factor = self.rnap.scaling_factor,
                                lb =0,
                                ub=1)

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
        k = self.rnap.ktrans / reaction.nucleotide_length

        # Scaled version:
        # rnap_constraint_expr = reaction.scaled_net - RMi # scaled
        rnap_constraint_expr = reaction.net - k * RMi.unscaled


        self.add_constraint(kind=SynthesisConstraint, hook=reaction,
                            expr=rnap_constraint_expr, ub=0)



    def init_rnap_variables(self, free_ratio):
        """
        Adds Free and Total ribosome variables to the models
        :return:
        """

        # Free ribosomes
        self._rnapf = self.add_variable(FreeRibosomes,
                                    self.rnap,
                                    lb=0,
                                    ub=1,
                                    scaling_factor = self.rnap.scaling_factor)

        # Add constraint on availability of free ribosomes
        # Rf and Rt are both concentrations
        expr = self.RNAPf - free_ratio * self.rnap.variable #scaled
        self.add_constraint(EnzymeRatio,
                            hook=self,
                            expr=expr,
                            id_='rnap',
                            lb=0,
                            ub=0)

    @property
    def RNAPf(self):
        return self._rnapf


    def add_ribosome(self, ribosome, free_ratio = 0.2):
        """
        Adds the ribosome used by the model.

        :param ribosome:
        :type ribosome: pytfa.me.Ribosome
        :return:
        """

        self.ribosome = ribosome

        self.add_enzymes(ribosome)

        # 0 -> We still need to add the virtual complexation of ribosomes:
        # it will be the same for all the translations of the model, so we can
        # call it from the ribosomal protein translation for example

        rprot_stoich = defaultdict(int)

        # Peptides
        for rprot_id in ribosome.composition:
            rprot_stoich[self.peptides.get_by_id(rprot_id)] -= 1


        complexation = ProteinComplexation(id='rib_complex',
                                           name='Ribosome complexation',
                                           target = ribosome,
                                           scaled=True)
        self.add_reactions([complexation])
        complexation.add_peptides(rprot_stoich)

        # Add it to a specific index
        self.complexation_reactions += [complexation]
        self.ribosome.init_variable()

        # Let us not forget the degradation
        self._add_enzyme_degradation(ribosome, scaled=True)

        self.init_ribosome_variables(free_ratio=free_ratio)

    def add_rrnas_to_rib_assembly(self):
        """
        Has to be done after the transcription reactions have been added
        :return:
        """
        # rRNA
        rrnas = []

        f = self.ribosome.kdeg * self.ribosome.scaling_factor
        # f = self.mu_max * self.ribosome.scaling_factor

        for the_rrna_id in self.ribosome.rrna_composition:
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


        self.ribosome.complexation.add_metabolites(
            {x:-1 for x in rrnas}, rescale=True)
            # {x:-1*f for x in rrnas}, rescale=False)
        self.rrnas += rrnas


    def init_ribosome_variables(self, free_ratio):
        """
        Adds Free and Total ribosome variables to the models
        :return:
        """

        # Free ribosomes
        self._Rf = self.add_variable(FreeRibosomes,
                                    self.ribosome,
                                    lb=0,
                                    ub=1,
                                    scaling_factor = self.ribosome.scaling_factor)

        # Add constraint on availability of free ribosomes
        # Rf and Rt are both concentrations
        expr = self.Rf - free_ratio * self.Rt #scaled
        self.add_constraint(EnzymeRatio,
                            hook=self,
                            expr=expr,
                            id_='rib',
                            lb=0,
                            ub=0)


    @property
    def Rt(self):
        return self.ribosome.variable
    @property
    def Rf(self):
        return self._Rf

    def _populate_ribosomes(self, free_ratio=0.2):
        """
        Once ribosomes have been assigned to the model, we still need to link
        them to the rest of the variables and constraints. This function creates
        the mass balance constraint on the ribosomes, as well as the total
        ribosome capacity constraint
        :return:
        """

        self.add_rrnas_to_rib_assembly()

        # v_complexation =   complexation.forward_variable  \
        #                  - complexation.reverse_variable

        # 1 -> Write the ribosome mass balance
        # Total amount of ribosome is in:
        # mass_balance_expr =   v_complexation            \
        #                     - self.ribosome.kdeg  * Rt  \
        #                     - self.mu             * Rt

        # Create the mass balance constraint
        self.add_mass_balance_constraint(self.ribosome.complexation,
                                         self.ribosome)

        # 2 -> Parametrize all the translation reactions with ribosomal vmax
        for trans_rxn in self.translation_reactions:
            self.apply_ribosomal_catalytic_constraint(trans_rxn)

        # 3 -> Add ribosomal capacity constraint
        self.regenerate_variables()

        # CATCH : This is summing ~1500+ variable objects, and for a reason
        # sympy does not like it. Let's cut it in smaller chunks and sum
        # afterwards
        # sum_RPs = sum(self.get_variables_of_type(RibosomeUsage))
        all_ribosome_usage = self.get_variables_of_type(RibosomeUsage)

        # sum_RPs = chunk_sum(all_ribosome_usage)
        sum_RPs = symbol_sum([x.unscaled for x in all_ribosome_usage])
        # For nondimensionalization
        ribo_usage = (sum_RPs + self.Rf.unscaled - self.ribosome.concentration) \
                     / self.ribosome.scaling_factor

        # Create the capacity constraint
        self.add_constraint(kind=TotalCapacity,
                            hook=self.ribosome,
                            expr=ribo_usage,
                            lb = 0,
                            ub = 0,
                            )

        # update variable and constraints attributes
        self.regenerate_constraints()
        self.regenerate_variables()


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

        RPi = self.add_variable(RibosomeUsage,
                                reaction.gene,
                                scaling_factor=self.ribosome.scaling_factor,
                                ub=1,
                                lb=0)

        # v_fwd - v_bwd <= kribo/length_aa [Ri]
        # v_fwd - v_bwd -  kribo/length_aa [Ri] <= 0

        # No scaling : Flux is in protein scale, and so is the ribosome concentration
        # v_max = self.ribosome.kribo \
        #         / reaction.aminoacid_length \
        #         * RPi
        # ribo_constraint_expr = fwd_variable - bwd_variable - v_max

        k = self.ribosome.kribo / reaction.aminoacid_length

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
        if hasattr(genes,'__iter__'):
            for g in genes:
                g._model = self
            self.genes += genes
        else:
            genes._model = self
            self.genes += [genes]


    #-------------------------------------------------------------------------#

    def sanitize_varnames(self):
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

        :param memo:
        :return:
        """

        return self.copy()

    def copy(self):

        from ..io.dict import model_from_dict, model_to_dict
        dictmodel = model_to_dict(self)
        new = model_from_dict(dictmodel)

        copy_solver_configuration(self, new)

        return new
