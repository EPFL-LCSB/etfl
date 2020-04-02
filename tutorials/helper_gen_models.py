#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#       E. Coli. Model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import logging

import pandas as pd

from cobra.flux_analysis.variability import flux_variability_analysis

from pytfa.io import        load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data
from pytfa.optim.relaxation import relax_dgo

from etfl.core import Enzyme, Ribosome, RNAPolymerase, ThermoMEModel, MEModel
from etfl.core.allocation import add_protein_mass_requirement, \
    add_rna_mass_requirement, add_dna_mass_requirement

from etfl.io.json import save_json_model
from pytfa.utils.logger import get_timestr

from etfl.data.ecoli import   get_model, get_thermo_data, get_coupling_dict, \
                        get_mrna_dict, get_rib, get_rnap, get_monomers_dict, \
                        get_nt_sequences, get_ratios, get_neidhardt_data, \
                        get_mrna_metrics, get_enz_metrics, \
                        remove_from_biomass_equation, get_ecoli_gen_stats, \
                        get_lloyd_coupling_dict, get_transporters_coupling, \
                        get_essentials, get_average_kcat, get_dna_polymerase

from etfl.optim.config import standard_solver_config
from etfl.analysis.summary import print_standard_sol

from optlang.exceptions import SolverError

from multiprocessing import Pool


# Run model gen in parallel ?
PARALLEL = False

data_dir = '../organism_data/info_ecoli'

# solver = 'optlang-cplex'
solver = 'optlang-gurobi'

# McCloskey2014 values
glc_uptake = 7.54
glc_uptake_std = 0.56
observed_growth_std = 0.02
observed_growth = 0.61

growth_reaction_id = 'BIOMASS_Ec_iJO1366_WT_53p95M'

def apply_bounds(model, bounds):
    for rid in bounds.index:
        rxn = model.reactions.get_by_id(rid)
        rxn.lower_bound = bounds.loc[rid].min()
        rxn.upper_bound = bounds.loc[rid].max()

# Free RNAP ratio should be around 0.7
# https://www.sciencedirect.com/science/article/pii/S0300908403001056?via%3Dihub
# Biochimie
# Volume 85, Issue 6, June 2003, Pages 597-609
# Free RNA polymerase and modeling global transcription in Escherichia coli
# H Bremer, P Dennis, M Ehrenberg

def create_model(has_thermo, has_expression, has_allocation, 
				kcat_mode='kmax',
                infer_missing_enz=False,
                additional_enz = None,
				free_rib_ratio=0.2,
				free_rnap_ratio=0.75,
				add_displacement = False,
				n_mu_bins = 128,
                name_suffix='',
                kcat_overrides=None):
    #------------------------------------------------------------
    # Initialisation
    #------------------------------------------------------------

    assert has_expression == True

    # this hack works because we are using the solver switch to update the var
    # names in the solver but really we should not do this
    # TODO: clean up model.sanitize_varnames
    vanilla_model = get_model('optlang-glpk')
    vanilla_model.reactions.EX_glc__D_e.lower_bound = -1 * glc_uptake - glc_uptake_std
    vanilla_model.reactions.EX_glc__D_e.upper_bound = -1 * glc_uptake + glc_uptake_std
    vanilla_model.objective = growth_reaction_id
    fba_sol = vanilla_model.slim_optimize()

    # vanilla_model.reactions.get_by_id(growth_reaction_id).lower_bound = observed_growth
    # fva = flux_variability_analysis(vanilla_model)
    # vanilla_model.reactions.get_by_id(growth_reaction_id).lower_bound = 0

    # original_bounds = pd.DataFrame.from_dict(
    #     {r.id:(r.lower_bound, r.upper_bound)
    #      for r in vanilla_model.reactions}, orient = 'index')
    # original_bounds.columns = ['lb','ub']


    mu_0 = fba_sol
    mu_range = [0, 3.5]
    n_mu_bins = n_mu_bins

    time_str = get_timestr()

    coupling_dict = get_coupling_dict(vanilla_model,
                                      mode = kcat_mode,
                                      atps_name = 'ATPS4rpp',
                                      infer_missing_enz = infer_missing_enz)


    if additional_enz is not None:
        additional_dict = get_transporters_coupling(model = vanilla_model,
                                                    additional_enz=additional_enz)

        additional_dict.update(coupling_dict)

        coupling_dict = additional_dict

    if kcat_overrides is not None:
        for rxn,enz_list in coupling_dict.items():
            for e in enz_list:
                if e.id in kcat_overrides:
                    prev_kcat = e.kcat_fwd
                    new_kcat = kcat_overrides[e.id]
                    e.kcat_fwd = new_kcat
                    e.kcat_bwd = new_kcat
                    print('Replaced kcat for {}: {} <-- {} s-1'
                          .format(e.id, prev_kcat/3600,new_kcat/3600))


    aa_dict, rna_nucleotides, rna_nucleotides_mp, dna_nucleotides = get_monomers_dict()
    essentials = get_essentials()

    # Initialize the model
    model_name = 'ETFL' if has_thermo else 'EFL'
    model_name = ('v'+model_name) if has_allocation else model_name
    model_name = (model_name + '_{}'.format(name_suffix)) if name_suffix else model_name
    model_name = (model_name+'_infer') if bool(infer_missing_enz) else model_name

    name = 'iJO1366_{}_{}_enz_{}_bins_{}.json'.format(
        model_name,
        len(coupling_dict),
        n_mu_bins,
        time_str)


    if has_thermo:

        thermo_data, lexicon, compartment_data = get_thermo_data()

        ecoli = ThermoMEModel(thermo_data,model = vanilla_model,
                              growth_reaction = growth_reaction_id,
                              mu_range = mu_range,
                              n_mu_bins = n_mu_bins,
                              name = name,
                              )
    else:
        ecoli = MEModel(model = vanilla_model,
                        growth_reaction = growth_reaction_id,
                        mu_range = mu_range,
                        n_mu_bins = n_mu_bins,
                        name = name,
                        )

    ecoli.name = name
    ecoli.logger.setLevel(logging.WARNING)
    ecoli.sloppy = True
    # apply_bounds(ecoli,fva)

    ecoli.solver = solver
    standard_solver_config(ecoli, verbose=False)

    if has_thermo:
        # Annotate the cobra_model
        annotate_from_lexicon(ecoli, lexicon)
        apply_compartment_data(ecoli, compartment_data)

        # TFA conversion
        ecoli.prepare()
        ecoli.convert(add_displacement = add_displacement)


    mrna_dict = get_mrna_dict(ecoli)
    nt_sequences = get_nt_sequences()
    rnap = get_rnap()
    # rnap.kcat_fwd *= 0.5
    rib = get_rib()

    # Remove nucleotides and amino acids from biomass reaction as they will be
    # taken into account by the expression

    remove_from_biomass_equation(model = ecoli,
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

    ecoli.add_nucleotide_sequences(nt_sequences)
    ecoli.add_essentials(  essentials=essentials,
                           aa_dict=aa_dict,
                           rna_nucleotides=rna_nucleotides,
                           rna_nucleotides_mp = rna_nucleotides_mp
                           )
    ecoli.add_mrnas(mrna_dict.values())
    ecoli.add_ribosome(rib, free_rib_ratio)
    ecoli.add_rnap(rnap, free_rnap_ratio)

    ecoli.build_expression()
    ecoli.add_enzymatic_coupling(coupling_dict)

    if has_allocation:

        nt_ratios, aa_ratios = get_ratios()
        chromosome_len, gc_ratio = get_ecoli_gen_stats()
        kdeg_mrna, mrna_length_avg  = get_mrna_metrics()
        kdeg_enz,  peptide_length_avg   = get_enz_metrics()
        neidhardt_mu, neidhardt_rrel, neidhardt_prel, neidhardt_drel = get_neidhardt_data()

        ecoli.add_dummies(nt_ratios=nt_ratios,
                          mrna_kdeg=kdeg_mrna,
                          mrna_length=mrna_length_avg,
                          aa_ratios=aa_ratios,
                          enzyme_kdeg=kdeg_enz,
                          peptide_length=peptide_length_avg)
        add_protein_mass_requirement(ecoli, neidhardt_mu, neidhardt_prel)
        add_rna_mass_requirement(ecoli, neidhardt_mu, neidhardt_rrel)
        add_dna_mass_requirement(ecoli, mu_values=neidhardt_mu,
                                       dna_rel=neidhardt_drel,
                                       gc_ratio=gc_ratio,
                                       chromosome_len=chromosome_len,
                                       dna_dict=dna_nucleotides)

        dna_pol = get_dna_polymerase()
        ecoli.add_enzymatic_coupling({'DNA_formation':[dna_pol,]})

    # Need to put after, because dummy has to be taken into account if used.
    ecoli.populate_expression()
    ecoli.add_trna_mass_balances()


    ecoli.print_info()
    ecoli.growth_reaction.lower_bound = observed_growth - 1*observed_growth_std

    need_relax = False

    ecoli.repair()

    try:
        ecoli.optimize()
    except (AttributeError, SolverError):
        need_relax = True

    if has_thermo and need_relax:
        # final_model, slack_model, relax_table = relax_dgo(ecoli)
        final_model, slack_model, relax_table = relax_dgo(ecoli, in_place = True)
    else:
        final_model = ecoli


    final_model.growth_reaction.lower_bound = 0
    # apply_bounds(ecoli, original_bounds)
    solution = final_model.optimize()
    print_standard_sol(final_model)

    filepath = 'models/{}'.format(final_model.name)
    save_json_model(final_model, filepath)

    final_model.logger.info('Build complete for model {}'.format(final_model.name))

    return final_model


if __name__ == '__main__':

    lloyd_enz = ['G1PPpp', 'GLCDpp', 'GLCabcpp', 'GLCptspp', 'GLCt2pp',
                 'LACZpp', 'LCTSt3ipp', 'LCTStpp','LACZ',
                 # 'ACt4pp','ACt2rpp',
                 'GALKr','GALt2pp','GALabcpp',
                 'ICL','MALS','PPCK','PPS','PGMT',
                 ]

    models = dict()

    # Models defined by Thermo - Expression - Neidhardt
    model_calls = dict()
    # model_calls[ 'EFL'  ] = {'has_expression':True,
    #                          'has_thermo':False,
    #                          'has_allocation':False,
    #                          'has_allocation':False,
    #                          'kcat_mode':'kmax',
    #                          'name_suffix':'v_0.11'}
    # model_calls[ 'ETFL' ] = {'has_expression':True,
    #                          'has_thermo':True,
    #                          'has_allocation':False,
    #                          'kcat_mode':'kmax',
    #                          'name_suffix':'v_0.11'}
    # model_calls['vEFL'  ] = {'has_expression':True,
    #                          'has_thermo':False,
    #                          'has_allocation':True,
    #                          'kcat_mode':'kmax',
    #                          'name_suffix':'v_0.11'}
    # model_calls['vETFL' ] = {'has_expression':True,
    #                          'has_thermo':True,
    #                          'has_allocation':True,
    #                          'kcat_mode':'kmax',
    #                          'name_suffix':'v_0.11'}
    # model_calls['vETFL65' ] = {'has_expression':True,
    #                            'has_thermo':True,
    #                            'has_allocation':True,
    #                            # 'kcat_mode':65*3600}
    #                            'kcat_mode':171.7*3600,
    #                            'name_suffix':'mean_kcat'}
    # model_calls['vETFL_infer' ]   = {   'has_expression':True,
    #                                     'has_thermo':True,
    #                                     'has_allocation':True,
    #                                     'kcat_mode':'kmax',
    #                                     'infer_missing_enz':True,
    #                                'name_suffix':'infer'}
    # model_calls['vETFL65_infer' ] = {   'has_expression':True,
    #                                     'has_thermo':True,
    #                                     'has_allocation':True,
    # #                                     # 'kcat_mode': 65*3600,
    #                                     'kcat_mode': 171.7*3600,
    #                                     'infer_missing_enz':True,
    #                                     'name_suffix':'infer_mean_kcat'}

    model_calls['vETFL_tp' ] = {'has_expression':True,
                                'has_thermo':True,
                                'has_allocation':True,
                                'kcat_mode':'kmax',
                                'name_suffix':'_tp_v_0.11',
                                'additional_enz':lloyd_enz}

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
