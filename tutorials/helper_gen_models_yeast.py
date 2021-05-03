#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#       Yeast. Model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import logging

import pandas as pd

from cobra.flux_analysis.variability import flux_variability_analysis

from pytfa.io import        load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data
from pytfa.optim.relaxation import relax_dgo

from etfl.core import Enzyme, Ribosome, RNAPolymerase, ThermoMEModel, MEModel

from etfl.io.json import save_json_model
from pytfa.utils.logger import get_timestr

from etfl.organisms.yeast import   get_model, get_thermo_data, get_coupling_dict, \
                        get_mrna_dict, get_rnap, get_monomers_dict, \
                        get_nt_sequences, get_ratios, get_mass_ratios, \
                        get_mrna_metrics, get_enz_metrics, \
                        remove_from_biomass_equation, get_yeast_gen_stats, \
                        get_essentials,  constrain_enzymes, get_rib, \
                        get_aa_sequences, update_copy_numbers, trna_charging_enz, \
                        get_rib_mit, get_rnap_mit, get_transcription_dict, \
                        get_translation_dict, remove_DNA, get_macromole_ratio, \
                        modify_GAM

from etfl.optim.config import standard_solver_config

from optlang.exceptions import SolverError

from multiprocessing import Pool
from etfl.debugging.debugging import relax_catalytic_constraints

from etfl.core.allocation import add_protein_mass_requirement, \
	    add_rna_mass_requirement, add_dna_mass_requirement, \
            add_lipid_mass_requirement, add_carbohydrate_mass_requirement,\
                add_ion_mass_requirement, fix_prot_ratio, fix_RNA_ratio, fix_DNA_ratio
from etfl.analysis.summary import print_standard_sol
from etfl.optim.constraints import RNAPAllocation

from etfl.organisms.yeast_utils import coupling_trna_enzymes_dict

# Costrain further the enzymes just like gecko?
frac_proteome = True

# Run model gen in parallel ?
PARALLEL = False # True

data_dir = '../organism_data/info_yeast/'

#solver='optlang-glpk'
#solver = 'optlang-cplex'
solver = 'optlang-gurobi'

# Martinez2014 values
glc_uptake = 15.2
glc_uptake_std = 0.1
observed_growth_std = 0.01
observed_growth = 0.43

growth_reaction_id = 'r_4041' #yeast 8 biomass pseudoreaction

def apply_bounds(model, bounds):
    for rid in bounds.index:
        rxn = model.reactions.get_by_id(rid)
        rxn.lower_bound = bounds.loc[rid].min()
        rxn.upper_bound = bounds.loc[rid].max()


def create_model(has_thermo, has_expression, var_allocation, 
				kcat_mode='kmax',
                infer_missing_enz=False,
                additional_enz = None,
				free_rib_ratio=0.15, # Bionumber:103025
				free_rnap_ratio=0.5, 
				add_displacement = False,
				n_mu_bins = 128):
    #------------------------------------------------------------
    # Initialisation
    #------------------------------------------------------------

    assert has_expression == True

    # this hack works because we are using the solver switch to update the var
    # names in the solver but really we should not do this
    # clean up model.sanitize_varnames
    vanilla_model = get_model('optlang-glpk')
    # D-glucose exchange
    vanilla_model.reactions.r_1714.lower_bound = -1 * glc_uptake - glc_uptake_std
    vanilla_model.reactions.r_1714.upper_bound = -1 * glc_uptake + glc_uptake_std
    
#    vanilla_model.reactions.r_4046.lower_bound = 0 # non-growth maintenance association
    
    vanilla_model.objective = growth_reaction_id
    fba_sol = vanilla_model.slim_optimize()
    
    
    # we need to load an unmodiifed FBa model and extract the biomass composition
    # this data is used to define the allocation constraints on the model for the basic case
    mass_ratios = get_macromole_ratio()
    mass_ratios['lipid'] = 0.07 # based on GECKO
    modify_GAM(vanilla_model, growth_reaction_id, prot_rxn_id='r_4047')

    mu_0 = fba_sol
    mu_range = [0, 1.5]
    n_mu_bins = n_mu_bins

    time_str = get_timestr()

    coupling_dict = get_coupling_dict(vanilla_model,
                                      mode = kcat_mode,
                                      atps_name = 'r_0226',
                                      infer_missing_enz = infer_missing_enz)

    # if additional_enz is not None:
        



    aa_dict, rna_nucleotides, rna_nucleotides_mp, dna_nucleotides = get_monomers_dict()
    essentials = get_essentials()

    # Initialize the model
    model_name = 'ETFL' if has_thermo else 'EFL'
    model_name = ('v'+model_name) if var_allocation else ('c'+model_name)
    model_name = (model_name+'_infer') if bool(infer_missing_enz) else model_name

    name = 'yeast8_{}_{}_enz_{}_bins_{}'.format(
        model_name,
        len(coupling_dict),
        n_mu_bins,
        time_str)

    if has_thermo:

        thermo_data, lexicon, compartment_data = get_thermo_data()

        yeast = ThermoMEModel(thermo_data,model = vanilla_model,
                              growth_reaction = growth_reaction_id,
                              mu_range = mu_range,
                              n_mu_bins = n_mu_bins,
                              name = name,
                              )
    else:
        yeast = MEModel(model = vanilla_model,
                        growth_reaction = growth_reaction_id,
                        mu_range = mu_range,
                        n_mu_bins = n_mu_bins,
                        name = name,
                        )
    yeast.name = name
    yeast.logger.setLevel(logging.WARNING)
    yeast.sloppy = False
    # apply_bounds(yeast,fva)

    yeast.solver = solver
    standard_solver_config(yeast, verbose=False)

    if has_thermo:
        # Annotate the cobra_model
        annotate_from_lexicon(yeast, lexicon)
        apply_compartment_data(yeast, compartment_data)

        # TFA conversion
        yeast.prepare()
        
        # if we need to specifically relax a thermo constraint
#        relax_thermo(yeast)
        yeast.convert(add_displacement = add_displacement)
        # removing delta G for membrane associated reactions
        rxnID_DG = [x.name.replace('DGo_','') for x in yeast.variables \
                        if 'DGo_' in x.name]
        for rxn_id in rxnID_DG:
            rxn = yeast.reactions.get_by_id(rxn_id)
            met_comps = [yeast.compartments[met.compartment]['name'] for met in rxn.metabolites]
            mem_asscd = False if len([comp for comp in met_comps \
                                      if 'membrane' in comp]) == 0 else True
            if mem_asscd:
                cstr_id = 'G_' + rxn_id
                yeast.remove_constraint(yeast._cons_dict[cstr_id])
                var_id = 'DG_' + rxn_id
                yeast.remove_variable(yeast._var_dict[var_id])
                var_id = 'DGo_' + rxn_id
                yeast.remove_variable(yeast._var_dict[var_id])


    mrna_dict = get_mrna_dict(yeast)
    nt_sequences = get_nt_sequences()
    aa_sequences = get_aa_sequences()
    rnap = get_rnap()
    rib_a, rib_b = get_rib()
    rnap_mit = get_rnap_mit()
    rib_mit = get_rib_mit()
    
    transcription_dict = get_transcription_dict()
    translation_dict = get_translation_dict()

    # Remove nucleotides and amino acids from biomass reaction as they will be
    # taken into account by the expression
    # in the case of yeast8 this function does nothing!

    remove_from_biomass_equation(model = yeast,
                                 nt_dict = rna_nucleotides,
                                 aa_dict = aa_dict,
                                 atp_id=essentials['atp'],
                                 adp_id=essentials['adp'],
                                 pi_id=essentials['pi'],
                                 h2o_id=essentials['h2o'],
                                 h_id=essentials['h'],
                                 )

    ##########################
    ##    MODEL CREATION    ##
    ##########################
            
    yeast.add_nucleotide_sequences(nt_sequences)
    yeast.add_peptide_sequences(aa_sequences)
    yeast.add_essentials(  essentials=essentials,
                           aa_dict=aa_dict,
                           rna_nucleotides=rna_nucleotides,
                           rna_nucleotides_mp = rna_nucleotides_mp
                           )
    yeast.add_mrnas(mrna_dict.values())
    yeast.add_ribosome(rib_a,free_ratio=free_rib_ratio)
    yeast.add_ribosome(rib_b,free_ratio=free_rib_ratio)
    yeast.add_rnap(rnap, free_ratio=free_rnap_ratio)
    yeast.add_ribosome(rib_mit,free_ratio=free_rib_ratio)
    yeast.add_rnap(rnap_mit, free_ratio=free_rnap_ratio)
    # should be after adding ribosome, since the rrna genes are replaced
    yeast.add_transcription_by(transcription_dict)
    yeast.add_translation_by(translation_dict)

    yeast.build_expression()
    yeast.add_enzymatic_coupling(coupling_dict)
    
    # Add enzymes for tRNA charging reactions
    enz_list, enz_dict = trna_charging_enz(yeast)
    yeast.add_enzymes(enz_list)
    
    
    # Dummy protein and mRNA should be added anyway
    nt_ratios, aa_ratios = get_ratios()
    chromosome_len, gc_ratio = get_yeast_gen_stats()
    kdeg_mrna, mrna_length_avg  = get_mrna_metrics()
    kdeg_enz,  peptide_length_avg   = get_enz_metrics()
    yeast.add_dummies(nt_ratios=nt_ratios,
                          mrna_kdeg=kdeg_mrna,
                          mrna_length=mrna_length_avg,
                          aa_ratios=aa_ratios,
                          enzyme_kdeg=kdeg_enz,
                          peptide_length=peptide_length_avg)
            
    if var_allocation:

        remove_DNA(yeast) 
        
        mu, rna, prot, dna, lipid, carbohydrate, ion = get_mass_ratios()

        add_protein_mass_requirement(yeast, mu, prot)
        add_rna_mass_requirement(yeast, mu, rna)
        add_dna_mass_requirement(yeast, mu_values=mu,
                                       dna_rel=dna,
                                       gc_ratio=gc_ratio,
                                       chromosome_len=chromosome_len,
                                       dna_dict=dna_nucleotides,
                                       ppi='s_0633_c')
        add_lipid_mass_requirement(yeast, lipid_mets = ['s_1096_c'],
                                   mass_ratios = mass_ratios,
                                   mu_values=mu,
                                   lipid_rel=lipid,
                                   lipid_rxn='r_2108')
        add_carbohydrate_mass_requirement(yeast, carbohydrate_mets = ['s_3718_c'],
                                   mass_ratios = mass_ratios,
                                   mu_values=mu,
                                   carbohydrate_rel=carbohydrate,
                                   carbohydrate_rxn='r_4048')
        add_ion_mass_requirement(yeast, ion_mets = ['s_4206_c'],
                                   mass_ratios = mass_ratios,
                                   mu_values=mu,
                                   ion_rel=ion,
                                   ion_rxn='r_4599')
    else:
        # constant allocation constraint to fix the sahere of RNA and protein
        fix_prot_ratio(yeast, mass_ratios)
        fix_RNA_ratio(yeast, mass_ratios)
        # this is added to have a DNA variable anyway which is needed for the RNAP allocation
        # DNA should not be removed from the biomass in this case
        fix_DNA_ratio(yeast, mass_ratios=mass_ratios, gc_ratio=gc_ratio, 
                      chromosome_len=chromosome_len)

        
    update_copy_numbers(yeast)   
    # Need to put after, because dummy has to be taken into account if used.
    yeast.populate_expression()
    yeast.add_trna_mass_balances()
    
    # Assign enzymes, which have been added, to the tRNA charging reactions
    new_coupling_dict = coupling_trna_enzymes_dict(yeast, enz_dict)
    yeast.add_enzymatic_coupling(new_coupling_dict)
    
    # after trna enzymes to account for them
    if frac_proteome:
        constrain_enzymes(yeast, mass_ratios)


    yeast.print_info()
    yeast.growth_reaction.lower_bound = observed_growth # - 1*observed_growth_std

    need_relax = False

    yeast.repair()

    try:
        yeast.optimize()
    except (AttributeError, SolverError):
        need_relax = True

    
    if has_thermo and need_relax:
        # final_model, slack_model, relax_table = relax_dgo(yeast)
        yeast.solver.configuration.timeout = 10*3600 # increasing time limit as it's more difficult
        final_model, slack_model, relax_table = relax_dgo(yeast, in_place = True)
    else:
        final_model = yeast


    final_model.growth_reaction.lower_bound = 0
    final_model.solver.problem.Params.FeasibilityTol = 1e-9
#    activator_states, final_model, relaxation=relax_catalytic_constraints(final_model, min_growth= 0.5)
    # apply_bounds(yeast, original_bounds)
    solution = final_model.optimize()
    flux_dict = {'  Glucose uptake    ' : 'r_1714',
                 '  Growth            ' : final_model.growth_reaction.id}
    print_standard_sol(final_model, flux_dict = flux_dict)

    filepath = 'models/{}'.format(final_model.name)
    save_json_model(final_model, filepath)

    final_model.logger.info('Build complete for model {}'.format(final_model.name))

    return final_model

if __name__ == '__main__':

    

    models = dict()

    # Models defined by Thermo - Expression - Neidhardt
    model_calls = dict()
    
    
    model_calls[ 'cEFL'  ] = {'has_expression':True,
                              'has_thermo':False,
                              'var_allocation':False,
                              'kcat_mode':'kmax'}
#    model_calls[ 'cETFL' ] = {'has_expression':True,
#                               'has_thermo':True,
#                               'var_allocation':False,
#                               'kcat_mode':'kmax'}
#    model_calls['vEFL'  ] = {'has_expression':True,
#                               'has_thermo':False,
#                               'var_allocation':True,
#                               'kcat_mode':'kmax'}
#    model_calls['vETFL' ] = {'has_expression':True,
#                              'has_thermo':True,
#                              'var_allocation':True,
#                              'kcat_mode':'kmax'}


    if not PARALLEL:
        for mid,mc in model_calls.items():
            models[mid] = create_model(**mc)
    else:
        pool = Pool()

        for mid,mc in model_calls.items():
            def this_callback(result, mid=mid):
                models[mid] = result
            pool.apply_async(create_model, [], mc, callback=this_callback)

        pool.close()
        pool.join()

    print('Completed')


    # Make thermo model
    # make_thermo_model()

    # Save FBA model
    # make_fba_model()
