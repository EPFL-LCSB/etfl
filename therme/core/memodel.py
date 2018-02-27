"""
.. module:: thermome
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Core for the ME-part

"""

import numpy as np
import optlang
import pandas as pd
import sympy
from cobra import Model, Reaction, Gene
from cobra.core import Solution, DictList
from collections import defaultdict


from ..utils.parsing import parse_gpr
from ..utils.utils import replace_by_enzymatic_reaction, replace_by_me_gene
from .genes import ExpressedGene
from .enzyme import Enzyme, Peptide
from .reactions import EnzymaticReaction, ProteinComplexation, \
    TranslationReaction, TranscriptionReaction
from .expression import build_trna_charging, \
    make_stoich_from_aa_sequence, make_stoich_from_nt_sequence
from ..optim.constraints import CatalyticConstraint, \
    ModelConstraint, EnzymeConstraint, MassBalance, TranslationConstraint, \
    GrowthCoupling, TotalCapacity, ExpressionCoupling, RibosomeRatio, \
    GrowthChoice, LinearizationConstraint
from ..optim.variables import ModelVariable, GrowthActivation, \
    GeneVariable, LinearizationVariable, RibosomeUsage, RNAPUsage, FreeRibosomes
from pytfa.core.model import LCSBModel
from pytfa.optim.reformulation import petersen_linearization
from pytfa.optim.utils import chunk_sum, symbol_sum
from pytfa.utils.logger import get_bistream_logger



class MEModel(LCSBModel, Model):
    def __init__(self, model=Model(), growth_reaction='',
                 mu = None, mu_error = 0,
                 mu_range = None, n_mu_bins = 1,
                 max_enzyme_concentration = 1000,
                 big_M = 1000,
                 scaling = 1000,
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

        name = 'ME2-' + model.id if model.id else 'ME2 model'

        LCSBModel.__init__(self, model, name)

        self.logger = get_bistream_logger('ME model' + str(self.name))
        self.parent = model
        if model is not None:
            self.sanitize_varnames()

        self.max_enzyme_concentration = max_enzyme_concentration
        self.big_M = big_M


        self._var_dict = dict()
        self._cons_dict = dict()

        self.logger.info('# ME Model initialized')

        self._growth_reaction_id = growth_reaction

        self._mu_error = mu_error
        self._mu_range = mu_range
        self._n_mu_bins = n_mu_bins
        self._mu_in = mu

        self._scaling = scaling

        if mu is not None and mu_error == 0:
            self._mu = mu
        elif mu is not None and mu_error > 0:
            self._mu = optlang.Variable(id = 'mu', name = 'mu')
            self._mu.lb = mu-mu_error
            self._mu.ub = mu+mu_error
            self.init_mu_variables()
        elif mu_range is not None:
            self._mu = optlang.Variable(id = 'mu', name = 'mu')
            self._mu.lb = mu_range[0]
            self._mu.ub = mu_range[1]
            self._n_mu_bins = n_mu_bins
            self.init_mu_variables()
        else:
            # message = """ You need to supply either mu, or mu_range.
            #             If you supply mu_error, it must be positive."""
            message = "Empty model initialized"
            # raise ValueError(message)
            self.logger.info(message)

        self.aa_dict = dict()
        self.nt_dict = dict()
        self.trna_dict = dict()

        self.enzymes = DictList()
        self.mrnas = DictList()
        self.peptides = DictList()
        self.transcription_reactions = DictList()
        self.translation_reactions = DictList()
        self.complexation_reactions = DictList()


    @property
    def mu(self):
        return self._mu

    @mu.setter
    def mu(self, val, epsilon = None):
        if epsilon is None:
            epsilon = self.solver.configuration.tolerances.feasibility

        self._mu.lb = val-epsilon
        self._mu.ub = val+epsilon

    def make_mu_bins(self):
        from numpy import linspace
        bounds = linspace(self.mu.lb, self.mu.ub, self._n_mu_bins+1)
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

        the_integer = sum([(2 ** i) * ga_i for i, ga_i in enumerate(ga)])

        binarized_mu = self.mu.lb + the_integer * (self.mu.ub - self.mu.lb) / N

        growth_coupling_expr = v_fwd - v_bwd - binarized_mu

        self.add_constraint(kind=GrowthCoupling,
                            hook=self.growth_reaction,
                            expr=growth_coupling_expr,
                            ub=bin_half_width,
                            lb=-1 * bin_half_width)

        # Update the variable indices
        self.regenerate_variables()
        self.regenerate_constraints()

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


    def build_expression(self, aa_dict, nt_dict,
                         rnap_genes,
                         rrna_genes,
                         rprot_genes,
                         atp='atp_c',
                         amp='amp_c',
                         gtp='gtp_c',
                         gdp='gdp_c',
                         ppi='ppi_c',
                         h2o='h2o_c',
                         h='h_c',
                         ):
        """
        Given a dictionnary from amino acids nucleotides to metabolite names,
        goes through the list of genes in the model that have sequence
        information to build transcription and traduction reactions

        :param aa_dict: A dictionnary of aminoacid letter to amicoacid met id
            Example :
            ```python
            aa_dict = {
                        'A':'ala__L_c',
                        'R':'arg__L_c',
                        ...
                    }
            ```
        :param nt_dict: A dictionnary of RNA nucleotide letter to nucleotide met id
            Example :
            ```python
            nt_dict = {
                        'A':'ade_c',
                        'U':'ura_c',
                        ...
                    }
            ```
        :param atp: atp metabolite id in the model
        :param amp: amp metabolite id in the model
        :param gtp: gtp metabolite id in the model
        :param gdp: gdp metabolite id in the model
        :param ppi: ppi metabolite id in the model
        :param h2o: h2o metabolite id in the model
        :param h: proton metabolite id in the model
        :return:
        """

        self.aa_dict = aa_dict
        self.nt_dict = nt_dict
        self.rnap_genes = rnap_genes
        self.rrna_genes = rrna_genes
        self.rprot_genes = rprot_genes

        self.trna_dict = build_trna_charging(self,aa_dict,atp,amp,ppi,h2o,h)

        # Check that the ribosomes have been added
        if self.ribosome is None:
            raise Exception(
                'A ribosome has to be added with the add_ribosome method')

        # Check that the RNAP has been added
        if self.rnap is None:
            raise Exception(
                'A RNA Polymerase has to be added with the add_rnap method')

        for gene in self.genes:
            if not isinstance(gene, ExpressedGene):
                continue

            # Build the transcription
            self._add_gene_transcription_reaction(gene)
            # Build the translation
            self._add_gene_translation_reaction(gene,gtp,gdp,h2o,h)


    def _add_gene_translation_reaction(self, gene,gtp,gdp,h2o,h):
        """

        :param gene: A gene of the model that has sequence data
        :type gene: therme.core.ExpressedGene
        :return:
        """

        rxn = TranslationReaction(
            id='translation_{}'.format(gene.id),
            name='Translation, {}'.format(gene.id),
            gene= gene,
            enzyme=self.ribosome)
        self.add_reactions([rxn])

        aa_stoichiometry = make_stoich_from_aa_sequence(gene.peptide,
                                                        self,
                                                        self.aa_dict,
                                                        self.trna_dict,
                                                        gtp,
                                                        gdp,
                                                        h2o,
                                                        h
                                                        )

        # Scale the stoichiometry
        aa_stoichiometry_scaled = {k:v/self._scaling \
                                   for k,v in aa_stoichiometry.items()}

        rxn.add_metabolites(aa_stoichiometry_scaled)

        free_peptide = Peptide(id = gene.id,
                               name = 'Peptide, {}'.format(gene.id))

        rxn.add_metabolites({free_peptide:1})

        # Add ribosome as necessary enzyme
        rxn.gene_reaction_rule = self.ribosome.id


        self.translation_reactions += [rxn]
        self.peptides += [free_peptide]

        
    def _add_gene_transcription_reaction(self, gene):
        """

        :param gene: A gene of the model that has sequence data
        :type gene: therme.core.ExpressedGene
        :return:
        """
        rxn = TranscriptionReaction(
            id='transcription_{}'.format(gene.id),
            name='Transcription, {}'.format(gene.id),
            gene= gene,
            enzyme=self.rnap)
        self.add_reactions([rxn])

        nt_stoichiometry = make_stoich_from_nt_sequence(gene.rna,
                                                        self,
                                                        self.nt_dict)

        # Scale the stoichiometry
        nt_stoichiometry_scaled = {k:v/self._scaling \
                                   for k,v in nt_stoichiometry.items()}

        rxn.add_metabolites(nt_stoichiometry_scaled)

        # Add rnap as necessary enzyme
        rxn.gene_reaction_rule = ' & '.join(self.rnap_genes)


        self.transcription_reactions += [rxn]


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
        for rid in self.coupling_dict:
            r = self.reactions.get_by_id(rid)

            # If the reaction is not compatible we do not try to apply constraints
            if not self.is_me_compatible(r):
                continue

            if isinstance(r, EnzymaticReaction) and r.id in coupling_dict:
                # This is a proper enzymatic reaction and we can directly apply
                # the constraint
                self.logger.debug('Applying catalytic constraint to {}'. \
                                 format(rid))
                r.add_enzymes(coupling_dict[r.id])
                self.apply_gpr_catalytic_constraint(r)
            elif not isinstance(r, EnzymaticReaction) and r.id in coupling_dict:
                # This reaction needs to be transformed to an EnzymaticReaction
                self.logger.debug('Transforming and applying catalytic constraint to {}'. \
                                 format(rid))
                #TODO : Add enzymatic_reaction dictlist ??
                enzyme = coupling_dict[r.id]
                enz_r = replace_by_enzymatic_reaction(self, r.id, enzyme)
                self.apply_gpr_catalytic_constraint(enz_r)
            else:
                self.logger.error('Could not find reaction {} in the coupling dictionnary'.format(r.id))

        # update variable and constraints attributes
        self.regenerate_constraints()
        self.regenerate_variables()

    def apply_gpr_catalytic_constraint(self, reaction):
        """
        Apply a catalytic constraint using a gene-enzymes reaction rule (GPR)

        :param reaction:
        :return:
        """

        complexation = self.add_complexation_from_gpr(reaction)

        v_max = dict()

        # Write v_max constraint
        fwd_variable = reaction.forward_variable
        bwd_variable = reaction.reverse_variable

        protein2isozyme_dict = self.match_enzymes_to_complexes(reaction.enzymes,
                                                               complexation)

        for e, (enz, comp) in enumerate(protein2isozyme_dict):
            # If the enzymes has the same kcat for both directions
            # v_fwd + v_bwd <= kcat [E]
            # v_fwd + v_bwd - kcat [E] <= 0

            v_max[e] =  ( enz.kcat / self._scaling )* enz.variable

            self.add_mass_balance_constraint(comp, enz)

            comp.enzyme = enz
            enz.complexation = comp

        enz_constraint_expr = fwd_variable + bwd_variable - sum(v_max.values())

        self.add_constraint(kind=CatalyticConstraint, hook=reaction,
                            expr=enz_constraint_expr, ub=0)


    def add_mass_balance_constraint(self, complexation, enzyme):
        """
        Adds a mass balance constraint of the type
        d[E]/dt = 0 <=> v_complexation - k_deg*[E] - μ*[E] = 0
        :param complexation:
        :param enzyme:
        :return:
        """
        # Add the mass_balance constraint
        v_complexation =    complexation.forward_variable \
                            -  complexation.reverse_variable

        # This is different if mu is a variable: we need to take care of the
        # bilinear constraint
        if isinstance(self.mu, optlang.Variable):
            # replace μ*E by z = sum(ga_i*μ_i*E), with ga_i binary variables
            # choosing between the mu_i
            z = self.linearize_me(enzyme)
            mass_balance_expr =     v_complexation \
                                    -   enzyme.kdeg * enzyme.variable \
                                    -   z

        else:
            # μ is fixed
            mass_balance_expr =     v_complexation \
                                    -   enzyme.kdeg * enzyme.variable \
                                    -   self.mu * enzyme.variable

        self.add_constraint(kind=MassBalance,
                            hook=enzyme,
                            expr=mass_balance_expr,
                            lb=0, ub=0)


    def is_me_compatible(self, reaction):
        # Test if the GPR is a proper one:
        this_gpr = reaction.gene_reaction_rule
        is_proper_gpr = bool(this_gpr) and this_gpr != '[]'

        sym_gpr = parse_gpr(this_gpr)

        ret = True

        if not is_proper_gpr:
            # Then we cannot constrain
            self.logger.warning('Improper GPR for {}'.format(reaction.id))
            ret = False

        # Check that all the genes participating in this gpr have a translation
        # reaction:
        is_translated = {x: x.name in self.translation_reactions
                         for x in sym_gpr.free_symbols}
        if not all(is_translated.values()):
            self.logger.warning(
                'Not all peptides in the GPR of {} are translated: {}'.format(
                    reaction.id, is_translated))
            ret = False

        return ret

    def linearize_me(self, enzyme):
        """
        Performs Petersen linearization on μ*E to stay an MILP problem

        :return:
        """

        E = enzyme.variable

        # ga_i is a binary variable for the binary expansion f the fraction on N
        # of the max growth rate
        ga_vars = self.get_variables_of_type(GrowthActivation)
        ga_vars = sorted(ga_vars, key=lambda x: x.ix)

        out_expr = self.mu.lb

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
                                          ub=self.max_enzyme_concentration)

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
                                    lb=cons.lb)

            out_expr += (2 ** i) * model_z_i * (self.mu.ub - self.mu.lb) \
                        / self.n_mu_bins

        return out_expr


    def add_complexation_from_gpr(self,reaction):
        """
        Logically parses the GPR to automatically find isozymes ( logical OR )
        and subunits ( logical AND ), and creates the necessary complexation
        reactions: 1 per isozyme, requiring the peptides of each subunit

        :param reaction:
        :type reaction: cobra.Reaction
        :return:
        """

        this_gpr = reaction.gene_reaction_rule

        sym_gpr = parse_gpr(this_gpr)

        if isinstance(sym_gpr, sympy.Symbol):
            # GPR of the type: '(gene0)'
            # Gene <=> Protein
            isozymes = [sym_gpr]
        elif isinstance(sym_gpr, sympy.And):
            # GPR of the type: '(gene0 & gene1)'
            # Subunits of one enzyme
            isozymes = [sym_gpr]
        elif isinstance(sym_gpr, sympy.Or):
            # GPR of the type: '(gene0 | gene1)', '((gene0 & gene1) | gene2)'
            # Two isozymes that are the arguments of the OR
            isozymes = sym_gpr.args

        complexation = []

        for e,this_isozyme in enumerate(isozymes):
            this_id = '{}_complex_{}'.format(reaction.id,e)
            this_name = '{} Complexation {}'.format(reaction.id,e)
            this_complexation = ProteinComplexation(id = this_id,
                                                    name = this_name)

            # TODO: Link stoichiometry information of the subunits
            if isinstance(this_isozyme, sympy.And):
                # this is a GPR with several subunits
                peptides = {self.peptides.get_by_id(x.name):-1 \
                            for x in this_isozyme.args}
            elif isinstance(this_isozyme, sympy.Symbol):
                # there is only one subunit
                peptides = {self.peptides.get_by_id(this_isozyme.name): -1}
            else:
                #The GPR has been incorrectly parsed
                self.logger.error('Incorrect parsing of {}'.format(isozymes))
                raise TypeError

            this_complexation.add_metabolites(peptides)

            complexation += [this_complexation]

        self.add_reactions(complexation)
        # Add it to a specific index
        self.complexation_reactions += complexation

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
        if len(enzyme_list) == 0:
            return None

        # First check whether the enzymes exist in the model
        enzyme_list = [x for x in enzyme_list if x.id not in self.enzymes]

        for enz in enzyme_list:
            enz._model = self
            enz.init_variable()

        for enz in enzyme_list:
            enz.variable.ub = self.big_M

        self.enzymes += enzyme_list


    def add_mrnas(self, mrna_list):
        """
        Adds an Enzyme object, or iterable of Enzyme objects, to the model
        :param enzyme_list:
        :type enzyme_list:Iterable(Enzyme) or Enzyme
        :return:
        """
        if not hasattr(mrna_list, '__iter__'):
            mrna_list = [mrna_list]
        if len(mrna_list) == 0:
            return None

        # First check whether the enzymes exist in the model
            mrna_list = [x for x in mrna_list if x.id not in self.mrnas]

        for mrna in mrna_list:
            mrna._model = self
            mrna.init_variable()

        for mrna in mrna_list:
            mrna.variable.ub = self.big_M

        self.mrnas += mrna_list

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
                             kcat = enzyme.kcat,
                             kdeg = enzyme.kdeg,
                             # name = enzyme.name + ' - Replicate {}'.format(e)
                             )
            new_enzymes.append(new_enz)

        self.add_enzymes(new_enzymes)
        return new_enzymes


    def populate_expression(self):
        self._populate_rnap()
        self._populate_ribosomes()

        # Add the coupling between mRNA availability and ribosome charging
        # The number of ribosomes assigned to a mRNA species is lower than
        # the number of such mRNA times the max number of ribosomes that can sit
        # on the mRNA:
        # [RPi] <= loadmax_i*[mRNAi]
        #
        # loadmax is : len(peptide_chain)/occupation(ribo)
        # "Their distance from one another along the mRNA is at least the size
        # of the physical footprint of a ribosome (≈20 nm, BNID 102320, 105000)
        # which is the length of about 60 base pairs (length of
        # nucleotide ≈0.3 nm, BNID 103777), equivalent to ≈20 aa."
        # "http://book.bionumbers.org/how-many-proteins-are-made-per-mrna-molecule/"
        #
        # hence:
        # [RPi] <= L_nt/Ribo_footprint * [mRNA]

        ribo_footprint_size = 60


        for the_mrna in self.mrnas:
            RPi = self._var_kinds[RibosomeUsage][the_mrna.id].variable
            mrna_var = the_mrna.variable

            polysome_size = len(the_mrna.gene.rna) / ribo_footprint_size
            expression_coupling = RPi - polysome_size * mrna_var
            self.add_constraint(kind=ExpressionCoupling,
                                hook = the_mrna,
                                expr=expression_coupling,
                                queue=True,
                                ub = 0)

        self.regenerate_variables()
        self.regenerate_constraints()

    def add_rnap(self, rnap):
        """
        Adds the RNA Polymerase used by the model.

        :param rnap:
        :type rnap: pytfa.me.Ribosome
        :return:
        """

        self.rnap = rnap

        # This is enough to add the translation, it will happen with the
        # ribosome pseudogene in add_gene_translation_reactions()
        rnap_pseudo_gene = Gene(rnap.id, name='RNA Polymerase pseudo-gene')
        self.add_genes(rnap_pseudo_gene)

        self.add_enzymes(rnap)

    def _populate_rnap(self):
        """
        Once ribosomes have been assigned to the model, we still need to link
        them to the rest of the variables and constraints. This function creates
        the mass balance constraint on the ribosomes, as well as the total
        ribosome capacity constraint
        :return:
        """
        # 0 -> We still need to add the virtual complexation of RNA Polymerases:
        peptide_stoich = defaultdict(int)
        for rprot_id in self.rnap_genes:
            peptide_stoich += self.peptides.get_by_id(rprot_id)

        complexation = ProteinComplexation(id='rnap_complex',
                                           name='RNA Polymerase complexation')
        complexation.add_metabolites(peptide_stoich)

        # v_complexation =   complexation.forward_variable  \
        #                  - complexation.reverse_variable

        # 1 -> Write the RNAP mass balance
        # Create the mass ba;ance constraint
        self.add_mass_balance_constraint(complexation, self.rnap)

        # 2 -> Parametrize all the transcription reactions with RNAP vmax
        for trans_rxn in self.translation_reactions:
            self.apply_rnap_catalytic_constraint(trans_rxn)

        # 3 -> Add RNAP capacity constraint
        self.regenerate_variables()

        all_rnap_usage = self.get_variables_of_type(RNAPUsage)
        sum_RMs = symbol_sum(all_rnap_usage)

        usage = sum_RMs - self.rnap.variable

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
        ribosome usage
        :param reaction: a TranscriptionReaction
        :type reaction: TranscriptionReaction
        :return:
        """

        # Check that we indeed have a transcription reaction
        assert(isinstance(reaction, TranscriptionReaction))

        RMi = self.add_variable(RNAPUsage, reaction)

        fwd_variable = reaction.forward_variable
        bwd_variable = reaction.reverse_variable

        # v_fwd - v_bwd <= kribo/length_aa [Ri]
        # v_fwd - v_bwd -  kribo/length_aa [Ri] <= 0

        v_max = self.rnap.ktrans / reaction.nucleotide_length * RMi

        ribo_constraint_expr = fwd_variable - bwd_variable - v_max


        self.add_constraint(kind=TranslationConstraint, hook=reaction,
                            expr=ribo_constraint_expr, ub=0)


    def add_ribosome(self, ribosome, free_ratio = 0.2):
        """
        Adds the ribosome used by the model.

        :param ribosome:
        :type ribosome: pytfa.me.Ribosome
        :return:
        """

        self.ribosome = ribosome

        self.add_enzymes(ribosome)

        self.init_ribosome_variables(free_ratio=free_ratio)


    def init_ribosome_variables(self, free_ratio=0.2):
        """
        Adds Free and Total ribosome variables to the models
        :return:
        """
        ## Add variables related to the ribosomes:
        # Total ribsomes
        # self.Rt = self.ribosome.variable # = v_syn^rib/(mu+k_deg^rib)
        # Is a property now
        # Free ribosomes
        self.Rf = self.add_variable(FreeRibosomes, self.ribosome)

        # Add constraint on availability of free ribosomes

        expr = self.Rf - free_ratio * self.Rt
        self.add_constraint(RibosomeRatio,
                             hook=self,
                             expr=expr,
                             id_='rib',
                             lb=0,
                             ub=0)


    @property
    def Rt(self):
        return self.ribosome.variable


    def _populate_ribosomes(self):
        """
        Once ribosomes have been assigned to the model, we still need to link
        them to the rest of the variables and constraints. This function creates
        the mass balance constraint on the ribosomes, as well as the total
        ribosome capacity constraint
        :return:
        """
        # 0 -> We still need to add the virtual complexation of ribosomes:
        # it will be the same for all the translations of the model, so we can
        # call it from the ribosomal protein translation for example

        rrna_stoich = defaultdict(int)
        for rrna_id in self.rrna_genes:
            rrna_stoich += self.mrnas.get_by_id(rrna_id)

        rprot_stoich = defaultdict(int)
        for rprot_id in self.rprot_genes:
            rprot_stoich += self.peptides.get_by_id(rprot_id)

        rib_stoich = rrna_stoich + rprot_stoich

        complexation = ProteinComplexation(id='rib_complex', name='Ribosome complexation')
        complexation.add_metabolites(rib_stoich)

        self.add_reactions(complexation)
        # Add it to a specific index
        self.complexation_reactions += complexation

        # v_complexation =   complexation.forward_variable  \
        #                  - complexation.reverse_variable

        # 1 -> Write the ribosome mass balance
        # Total amount of ribosome is in:
        # mass_balance_expr =   v_complexation            \
        #                     - self.ribosome.kdeg  * Rt  \
        #                     - self.mu             * Rt

        # Create the mass balance constraint
        self.add_mass_balance_constraint(complexation, self.ribosome)

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
        sum_RPs = symbol_sum(all_ribosome_usage)

        ribo_usage = sum_RPs + self.Rf - self.Rt
        # ribo_usage = sum_RPs + Rf - Rt
        # ribo_usage = sum_RPs + Rf - Rt

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

        RPi = self.add_variable(RibosomeUsage, reaction)

        fwd_variable = reaction.forward_variable
        bwd_variable = reaction.reverse_variable

        # v_fwd - v_bwd <= kribo/length_aa [Ri]
        # v_fwd - v_bwd -  kribo/length_aa [Ri] <= 0

        v_max = self.ribosome.kribo / reaction.aminoacid_length * RPi

        ribo_constraint_expr = fwd_variable - bwd_variable - v_max


        self.add_constraint(kind=TranslationConstraint, hook=reaction,
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
        n_enzymatic_reactions   = len([x for x in self.reactions   \
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
