#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 20:10:30 2020

@author: omid
"""
from .memodel import MEModel, RNAP_FOOTPRINT_SIZE, RIBO_FOOTPRINT_SIZE
from .allocation import fix_RNA_ratio, fix_prot_ratio, fix_vector_ratio, \
    MRNA_WEIGHT_CONS_ID, PROT_WEIGHT_CONS_ID, \
    DNA_WEIGHT_CONS_ID, MRNA_WEIGHT_VAR_ID, PROT_WEIGHT_VAR_ID, \
    DNA_WEIGHT_VAR_ID, DNA_FORMATION_RXN_ID, PROT_CONSTANT_CONS_ID, \
        DNA_CONSTANT_CONS_ID, RNA_CONSTANT_CONS_ID, ENZYME_CONSTANT_CONS_ID, \
    define_prot_weight_constraint, define_mrna_weight_constraint, \
    define_dna_weight_constraint, get_dna_synthesis_mets

from .genes import ExpressedGene, CodingGene  
from ..optim.constraints import tRNAMassBalance, InterpolationConstraint, \
    TotalCapacity, EnzymeRatio, ConstantAllocation, MinimalAllocation, MinimalCoupling
from ..optim.variables import InterpolationVariable, EnzymeVariable, RibosomeUsage, RNAPUsage

from pytfa.optim.utils import symbol_sum

DNA_VECTOR_ID = 'DNA vector'

class RecombModel(MEModel):

    def __init__(self, me_model, inplace = True):
        """
        Hack to subclass instantiated object:
        https://stackoverflow.com/questions/33463232/subclassing-in-python-of-instantiated-superclass

        :param me_model:
        :param inplace:
        """
        if not inplace:
            new = me_model.copy()
            self.__dict__ = new.__dict__
        else:
            # new = me_model
            self.__dict__ = me_model.__dict__

        # self._me_model = new
        self._has_transcription_changed = False
        self._has_translation_changed = False
        self.vectors = dict()
        self.vector_copy_number = dict()
        self._biomass_metabolites = dict()

    # def __getattr__(self,attr):
    #     """
    #     Hack to subclass instantiated object:
    #     https://stackoverflow.com/questions/33463232/subclassing-in-python-of-instantiated-superclass
    #
    #     :param attr:
    #     :return:
    #     """
    #     return getattr(self._me_model,attr)
    
    def add_vector(self, vector,
                   copy_number=1):
        '''
        

        Parameters
        ----------
        vector 
        copy_number : int, optional
            Vector copy number. The default is 1.

        Returns
        -------
        None.

        '''

        self.vectors[vector.id] = vector
        self.vector_copy_number[vector.id] = copy_number

        # Add the ribosome and RNAP of the vector if they exist
        if vector.rnap is not None:
            self.add_vector_RNAP(vector.rnap)
        if vector.ribosome is not None:
            self.add_vector_ribosome(vector.ribosome)

        # Add vector genes
        self.add_genes(vector.genes)
        for g in vector.genes: # This adds the peptides
            if isinstance(g, CodingGene):
                self._make_peptide_from_gene(g.id)
                vector.peptides += [self.peptides.get_by_id(g.id)]
            
        # Edit the gene copy number by the number of plasmids.
        for g in vector.genes:
            self.genes.get_by_id(g.id).copy_number = g.copy_number * copy_number

        self.express_genes(vector.genes)

        # Add the mRNAs
        self.add_mrnas(vector.mrna_dict.values())
        self._push_queue()

        # Constraint the polysomes and RNAP allocation
        for g in vector.genes:
            m = self.mrnas.get_by_id(g.id)
            self.add_mrna_mass_balance(m)

            the_transcription = self.get_transcription(m.id)
            self.apply_rnap_catalytic_constraint(the_transcription,
                                                 queue=False)

            if isinstance(g, CodingGene):
                the_translation = self.get_translation(m.id)
                self.apply_ribosomal_catalytic_constraint(the_translation)

            self.regenerate_variables()
            self.regenerate_constraints()

            self._constrain_polysome(m)
            self._constrain_polymerase(g)

        # Adding reaction-enzyme coupling
        # First add enzymes seperately, so if a reaction is not coupled to any reaction, it is still added.
        self.add_enzymes(vector.proteins) 
        self.add_reactions(vector.reactions)
        self.add_enzymatic_coupling(vector.coupling_dict)

        # recompute mRNA-dependent constraints
        self.recompute_trna_balances()
        self.recompute_transcription(vector.default_rnap, vector.genes)
        self.recompute_translation(vector.default_rib, vector.genes)

        # This needs to be after adding genes and mRNAs
        self.reallocation(vector)
        self._push_queue()                
    
    @property
    def biomass_metabolites(self):
        return self._biomass_metabolites
    
    @biomass_metabolites.setter
    def biomass_metabolites(self, met_dict):
        self._biomass_metabolites = met_dict
            
    def recompute_trna_balances(self):

        #1. Remove tRNA balances
        trna_mb_cons = self.get_constraints_of_type(tRNAMassBalance)
        for the_cons in trna_mb_cons:
            self.remove_constraint(the_cons)

        #2. Remake tRNA balances
        self.add_trna_mass_balances()

    def add_vector_RNAP(self, rnap):
        """
        Adds the vector's RNAP to the RNAP pool of the cell

        :param rnap:
        :return:
        """
        #TODO: Clean this up
        try:
            free_rnap_ratio = rnap.free_ratio
        except AttributeError:
            # Get it from the current free ratio constraint
            rnap_id = list(self.rnap.keys())[0]
            cons = self.get_constraints_of_type(EnzymeRatio).get_by_id(rnap_id)
            rnap_tot_var = self.rnap[rnap_id].variable
            # The Enzyme Ratio constraint looks like
            # [E_free] - ρ*[E_tot] = 0,
            # which is equivalent to
            # [E_free] = ρ*[E_tot]
            # So the linear coeff will be negative
            free_rnap_ratio = abs( 
                cons.constraint.get_linear_coefficients([rnap_tot_var])[rnap_tot_var]
                                    )

        # This adds the RNAP as an enzyme, and also enforces its free ratio
        self.add_rnap(rnap, free_ratio=free_rnap_ratio)

    def add_vector_ribosome(self, ribosome, free_rib_ratio=0):
        """
        Adds the vector's ribosome to the ribosome pool of the cell

        :param rnap:
        :return:
        """
        self.add_ribosome(ribosome, free_ratio=free_rib_ratio)
        
    def reallocation(self, vector, ppi='ppi_c'):
        """
        Add an allocation constraint for the vector DNA. Recompute RNA and protein
        allocation (to add vector mRNA and peptide to the LHS). The RHS of all
        allocations (including the stoichiometric coefficient of biomass components)
        is updated.
        :return:
        """
        interpolation_constraints = self.get_constraints_of_type(InterpolationConstraint)
        interpolation_variables = self.get_variables_of_type(InterpolationVariable)
        if not interpolation_variables: # constant allocation
            
            #0. update the implicit allocations with stoichiometric coefficients
            # rel_mass_ratio = length_of_vector * copy_number / length_of _DNA
            rel_mass_ratio = vector.len * self.vector_copy_number[vector.id] / self.dna.len
            rel_vector_ratio = rel_mass_ratio * self.biomass_composition['DNA']
            self._recompute_composition_dict(rel_vector_ratio)
            
            
            #1. update the explicit allocation constraints
            self._recompute_allocation_cstr(vector)
    
    
            #TODO: update it for the case that the vector is RNA vector
            #2. Add a constant allocation constraint for vector DNA
            fix_vector_ratio(self, vector, 
                             dna_ratio=self.biomass_composition[DNA_VECTOR_ID],
                             ppi=ppi) # This also adds DNA species for the vector
            
        
        else: # variable allocation
            raise NotImplementedError()
            #2. Remove previous allocation constraints (only for RNA and protein)
            # # mRNA
            # mrna_weight_def_cons = interpolation_constraints.get_by_id(MRNA_WEIGHT_CONS_ID)
            # mrna_weight_var = interpolation_variables.get_by_id(MRNA_WEIGHT_VAR_ID)
    
            # # Proteins
            # prot_weight_def_cons = interpolation_constraints.get_by_id(PROT_WEIGHT_CONS_ID)
            # prot_weight_var = interpolation_variables.get_by_id(PROT_WEIGHT_VAR_ID)
    
            # #DNA
            # # dna_weight_def_cons  = interpolation_constraints.get_by_id(DNA_WEIGHT_CONS_ID)
            # # dna_weight_var = interpolation_variables.get_by_id(DNA_WEIGHT_VAR_ID)
    
            # for the_cons in [mrna_weight_def_cons, prot_weight_def_cons]:
            #     self.remove_constraint(the_cons)
    
            # #3. Apply new allocation constraints
            # define_prot_weight_constraint(self,prot_weight_var)
            # define_mrna_weight_constraint(self,mrna_weight_var)
            
    def _recompute_composition_dict(self, vector_ratio):
        '''
        A function to find new dict of biomass composition after adding vector
        with a relative ratio (vector_ratio).
    
        Parameters
        ----------

        vector_ratio : TYPE
            DESCRIPTION.
    
    
        '''
        new_composition_dict = dict()
        try:
            tot_mass = self.biomass_composition['total mass']
        except KeyError:
            tot_mass = 1 # the default weight of biomass is 1 gDW
        
        mod_factor = tot_mass/(tot_mass + vector_ratio)
        for component, ratio in self.biomass_composition.items():
            new_composition_dict[component] = mod_factor * ratio
        
        new_composition_dict[DNA_VECTOR_ID] = vector_ratio * mod_factor
        # changing the stiochiometric coefficients of DNA, lipid, etc.
        self._stoich_coeff_from_mass_ratio(new_composition_dict)
        
        self.biomass_composition = new_composition_dict
        
        return 
    
    def _stoich_coeff_from_mass_ratio(self, new_alloc_data,
                                  lumped_comp = True):
        '''
        This is a function to find new stoichiometric coefficients for biomass reaction
        based on experimental data for mass ratios.
        
        inputs:
            new_alloc_data: a dict or DataFrame of different biomass components with
                their experimental mass ratios
            macromole_ratio: dict of mass ratios of different biomass components in
                GEM.
            lumped_comp: indicates if biomass metabolites are lumped into pseudometabolites
            met_macromol_dict: a dictionary to relate metabolite ids with each type
                of macromolecules, e.g. protein, RNA, etc. Keys are met ids and values
                are macromolecule names compatible with tags in alloc_data & macromole_ratio
            
        ouputs:
            new model with modified stoichiometric coefficients in its growth reaction.
        '''
        # RNA and protein are already removed from the biomass equation
        removed = ['RNA', 'protein']
        # finding the original stoichiometric coefficients for each metabolite
        biomass_rxn = self.growth_reaction
        met_macromol_dict = {k:v for k,v in self.biomass_metabolites.items() \
                             if v not in removed}
        macromole_ratio = {k:v for k,v in self.biomass_composition.items() \
                           if k not in removed}
        
        if lumped_comp:
            # nothing to do :)
            org_coeffs = {x : biomass_rxn.get_coefficient(x) for x \
                          in met_macromol_dict.keys()}
        else:
            # first, we should lump different metabolites to relate them to macromolecules
            raise NotImplementedError()
            
        new_coeffs = {k : v * new_alloc_data[met_macromol_dict[k]] /\
                      macromole_ratio[met_macromol_dict[k]] for k,v in org_coeffs.items()}
        
        # the new stoichiometric coefficient is added to the previous one
        # to avoid redundancy I should first subtract the old coefficient
        change_coeffs = {k : v - org_coeffs[k] for k,v in new_coeffs.items()}
        biomass_rxn.add_metabolites(change_coeffs)
        
        return
    
    def _recompute_allocation_cstr(self, vector):
        '''
        

        Returns
        -------
        None.

        '''
        comp_dict = self.biomass_composition
        
        allocation_cstr = self.get_constraints_of_type(ConstantAllocation)
        cstr_prot = allocation_cstr.get_by_id(PROT_CONSTANT_CONS_ID)
        cstr_rna = allocation_cstr.get_by_id(RNA_CONSTANT_CONS_ID)
        cstr_dna = allocation_cstr.get_by_id(DNA_CONSTANT_CONS_ID)
        cstr_enz = allocation_cstr.get_by_id(ENZYME_CONSTANT_CONS_ID)
        
        phi = cstr_enz.constraint.ub / cstr_prot.constraint.ub
        
        # DNA constraint expression does not have to be changed
        cstr_dna.constraint.lb = comp_dict['DNA']
        cstr_dna.constraint.ub = comp_dict['DNA']
        
        # removing enzyme, protein, RNA allocation constraint to be replaced by the new
        for the_cons in [cstr_enz, cstr_prot, cstr_rna]:
            self.remove_constraint(the_cons)
        
        fix_prot_ratio(self, comp_dict['protein'])        
        fix_RNA_ratio(self, comp_dict['RNA'])
        
        # A specific allocation constraint for enzymes
        # We assume that the additional protein takes share from both enzymes and dummy
        # This share is determined based on phi: 
        # phi * (concentraion of vector proteins) from the enzymes &
        # (1-phi) * (concentraion of vector proteins) from the dummy
        enz_vars = self.get_variables_of_type(EnzymeVariable)
        
        vector_enz = ['EZ_{}'.format(x.id) for x in vector.proteins]
        vec_expression = symbol_sum([x for x in enz_vars if x.name in vector_enz])
        
        dummy_var = enz_vars.get_by_id('dummy_enzyme')
        
        self.add_constraint(kind = ConstantAllocation, 
                                 hook = self, 
                                 expr = dummy_var + (1-phi) * vec_expression,
                                 id_ = ENZYME_CONSTANT_CONS_ID,
                                 lb = (1-phi) * comp_dict['protein']
                                 )
        
        return
    
    def recompute_translation(self, default_rib, genes):
        """

        :return:
        """
        if isinstance(default_rib, list) and len(default_rib) == 1:
            default_rib = default_rib[0]
        
        if isinstance(default_rib, str): # ribosome id
            rib_id = default_rib
        elif not hasattr(default_rib, '__iter__'): # ribosome object
            rib_id = default_rib.id
        elif isinstance(default_rib, dict):
            rib_id = None
        else:
            raise ValueError('The defualt ribosome of the vector is ambiguous')
            
        if not hasattr(genes, '__iter__'):
            genes = [genes] 
            
        gene_ids = [g.id for g in genes if isinstance(g, CodingGene)] # RibosomeUsage is only defined for CodingGene
        rib_usage = self.get_variables_of_type(RibosomeUsage)
        cstrs = self.get_constraints_of_type(TotalCapacity)
        
        if rib_id:
            add_usage = [rib_usage.get_by_id(id_) for id_ in gene_ids]
            change_expr = symbol_sum([x.variable for x in add_usage])
            rib_variable = self.get_variables_of_type(EnzymeVariable).get_by_id(rib_id)
            for the_cstr in cstrs:
                if rib_variable.variable in the_cstr.constraint.variables: # to find TotalCapacity constraints associated with this ribosome
                    old_expr = the_cstr.constraint.expression
                    new_expr = old_expr + change_expr
                    the_cstr.change_expr(new_expr)
            
        else:
            raise NotImplementedError()

    def recompute_transcription(self, default_rnap, genes):
        """

        :return:
        """

        if isinstance(default_rnap, list) and len(default_rnap) == 1:
            default_rnap = default_rnap[0]
        
        if isinstance(default_rnap, str): # ribosome id
            rnap_id = default_rnap
        elif not hasattr(default_rnap, '__iter__'): # ribosome object
            rnap_id = default_rnap.id
        elif isinstance(default_rnap, dict):
            rnap_id = None
        else:
            raise ValueError('The defualt RNAP of the vector is ambiguous')
            
        if not hasattr(genes, '__iter__'):
            genes = [genes] 
            
        gene_ids = [g.id for g in genes]
        rnap_usage = self.get_variables_of_type(RNAPUsage)
        cstrs = self.get_constraints_of_type(TotalCapacity)
        
        if rnap_id:
            add_usage = [rnap_usage.get_by_id(id_) for id_ in gene_ids]
            change_expr = symbol_sum([x.variable for x in add_usage])
            rnap_variable = self.get_variables_of_type(EnzymeVariable).get_by_id(rnap_id)
            for the_cstr in cstrs:
                if rnap_variable.variable in the_cstr.constraint.variables: # to find TotalCapacity constraints associated with this rnap
                    old_expr = the_cstr.constraint.expression
                    new_expr = old_expr + change_expr
                    the_cstr.change_expr(new_expr)
            
        else:
            raise NotImplementedError()
            
    def change_vector_protein_alloc(self, new_phi_p, vector_id):
        '''
        A function to change the fraction of total protein that can be allocated to the plasmid protein

        Parameters
        ----------
        new_phi_p : new value for the fraction of proteome that can be allocated to plasmid protein
        vector_id: the ID of the vector

        Returns
        -------
        None.

        '''
        # take the allocation constraint
        for cstr in self.get_constraints_of_type(ConstantAllocation):
            if ENZYME_CONSTANT_CONS_ID in cstr.name:
                the_allocation = cstr
                
        enz_vars = self.get_variables_of_type(EnzymeVariable)
        vector = self.vectors[vector_id]
        
        vector_enz = ['EZ_{}'.format(x.id) for x in vector.proteins]
        vec_expression = symbol_sum([x for x in enz_vars if x.name in vector_enz])
        
        dummy_var = enz_vars.get_by_id('dummy_enzyme')
        
        new_expr = dummy_var + (1-new_phi_p) * vec_expression
        the_allocation.change_expr(new_expr)
        # the lower bound should not be changed, as this is phi and not phi_plasmid
        # the_allocation.constraint.lb = (1-new_phi_p) * self.biomass_composition['protein']
        
        return
    
    def change_tcpt_basal_activity(self, new_basal_act, vector_id):
        '''
        A function to update constraints after changing the trancriptional minimal activity

        Parameters
        ----------
        new_basal_act : TYPE
            DESCRIPTION.
        vector_id : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        vector = self.vectors[vector_id]
        genes = vector.genes
        min_activity_cstr = self.get_constraints_of_type(MinimalAllocation)
        rnap_usage_vars = self.get_variables_of_type(RNAPUsage)
        
        for gene in genes:
            try:
                the_cstr = min_activity_cstr.get_by_id(gene.id)
                RNAPi_hat = rnap_usage_vars.get_by_id(gene.id)
                n_loci = gene.copy_number
                scaling_factor = self.dna.scaling_factor / RNAPi_hat.scaling_factor
                loadmax = len(gene.sequence) / RNAP_FOOTPRINT_SIZE
                new_expr = RNAPi_hat \
                     - new_basal_act * loadmax * n_loci * scaling_factor * \
                         self.dna.scaled_concentration
                the_cstr.change_expr(new_expr)
                gene.min_tcpt_activity = new_basal_act
            except KeyError:
                pass
            
    def change_tnsl_basal_activity(self, new_basal_act, vector_id):
        '''
        A function to update constraints after changing the translational minimal activity

        Parameters
        ----------
        new_basal_act : TYPE
            DESCRIPTION.
        vector_id : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        vector = self.vectors[vector_id]
        genes = vector.genes
        min_activity_cstr = self.get_constraints_of_type(MinimalCoupling)
        rib_usage_vars = self.get_variables_of_type(RibosomeUsage)
        
        for gene in genes:
            try: # only coding genes have this constraint
                the_cstr = min_activity_cstr.get_by_id(gene.id)
                the_mrna = self.mrnas.get_by_id(gene.id)
                RPi_hat = rib_usage_vars.get_by_id(the_mrna.id)
                mrna_hat = the_mrna.scaled_concentration
                scaling_factor = the_mrna.scaling_factor / RPi_hat.scaling_factor
                polysome_size = len(the_mrna.gene.rna) / RIBO_FOOTPRINT_SIZE
                new_expr = RPi_hat \
                              - new_basal_act * polysome_size * scaling_factor * mrna_hat
                the_cstr.change_expr(new_expr)
                gene.min_tnsl_activity = new_basal_act
            except KeyError:
                pass