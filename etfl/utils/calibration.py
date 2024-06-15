#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A script for all kinds of calibration

"""

from ..optim.constraints import TotalCapacity
from ..optim.variables import RNAPUsage, RibosomeUsage, mRNAVariable

def calib_plasmid_activity(model, rnap_frac, rib_frac, product_id):
    '''
    

    Parameters
    ----------
    model : RecombMeModel, without activity, i.e. omega = 0
    rnap_frac : dict
        The fraction of rnap that should be allocated to plasmid.
    rib_frac : dict
        The fraction of rib that should be allocated to plasmid.
    product_id : str
        The ID for the protein that is produced by the plasmid.

    Returns
    -------
    omega_tcp : dict
        DESCRIPTION.
    omega_tnl : dict
        DESCRIPTION.

    '''
    
    omega_tcp = dict()
    omega_tnl = dict()
    
    rnap_footprint_size = 40
    ribo_footprint_size = 60
    
    rnap_usage_id = 'RM_' + product_id
    ribo_usage_id = 'RP_' + product_id
    
    for copy_number, f_rnap in rnap_frac.items():
        try:
            f_ribo = rib_frac[copy_number]
        except KeyError:
            raise Exception('The copy number values for rnap and ribosome are not compatible')
            
        # adding the constraint f_rnap*Ez_rnap - RM_p_beta_lac = 0
        rnap_cstr = model.add_constraint(kind=TotalCapacity,
                                    hook=model,
                                    id_ = 'rnap_plasmid',
                                    expr= f_rnap * model.variables.EZ_rnap - \
                                        model.variables.get(rnap_usage_id),
                                    lb = 0,
                                    ub = 0,
                                    )
        # adding the constraint f_rib*Ez_rib - RP_p_beta_lac = 0
        ribo_cstr = model.add_constraint(kind=TotalCapacity,
                                    hook=model,
                                    id_ = 'rib_plasmid',
                                    expr= f_ribo * model.variables.EZ_rib - \
                                        model.variables.get(ribo_usage_id),
                                    lb = 0,
                                    ub = 0,
                                    )
        # solving the problem
        model.slim_optimize()
        
        # back-calculating omega
        gene = model.genes.get_by_id(product_id)
        polysome_size = len(gene.sequence) / rnap_footprint_size
        n_loci = copy_number
        scaled_conc = model.variables.DN_DNA.primal
        rnap_usage = model.get_variables_of_type(RNAPUsage).get_by_id(gene.id)
        scaling_factor = model.dna.scaling_factor / rnap_usage.scaling_factor
        omega_tcp[copy_number] = \
            rnap_usage.variable.primal/(scaling_factor*scaled_conc*polysome_size*n_loci)
        ##    
        mrna = model.mrnas.get_by_id(product_id)
        polysome_size = len(mrna.gene.rna) / ribo_footprint_size
        scaled_conc = model.get_variables_of_type(mRNAVariable).get_by_id(mrna.id)
        rib_usage = model.get_variables_of_type(RibosomeUsage).get_by_id(mrna.id)
        scaling_factor = mrna.scaling_factor / rib_usage.scaling_factor
        omega_tnl[copy_number] = \
            rib_usage.variable.primal/(scaling_factor*scaled_conc.variable.primal*polysome_size)
        
        # removing the constraints for the next round
        model.remove_constraint(rnap_cstr)
        model.remove_constraint(ribo_cstr)
    
    return omega_tcp, omega_tnl