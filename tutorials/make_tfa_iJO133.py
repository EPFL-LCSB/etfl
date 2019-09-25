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

from etfl.io.json import save_json_model
from pytfa.utils.logger import get_timestr

from etfl.data.ecoli import   get_model, get_thermo_data, get_coupling_dict, \
                        get_mrna_dict, get_rib, get_rnap, get_monomers_dict, \
                        get_nt_sequences, get_ratios, get_neidhardt_data, \
                        get_mrna_metrics, get_enz_metrics, \
                        remove_from_biomass_equation, get_ecoli_gen_stats, \
                        get_lloyd_coupling_dict, get_transporters_coupling, \
                        get_essentials, get_average_kcat

from etfl.optim.config import standard_solver_config

from optlang.exceptions import SolverError

from multiprocessing import Pool
from pytfa.thermo.tmodel import ThermoModel


# solver = 'optlang-cplex'
solver = 'optlang-gurobi'

# McCloskey2014 values
glc_uptake = 7.54
glc_uptake_std = 0.56
observed_growth_std = 0.02
observed_growth = 0.61

growth_reaction_id = 'BIOMASS_Ec_iJO1366_WT_53p95M'


def create_tfa_model(add_displacement = False):
    #------------------------------------------------------------
    # Initialisation
    #------------------------------------------------------------

    time_str = get_timestr()
    name = 'iJO1366_TFA_{}.json'.format(time_str)

    # this hack works because we are using the solver switch to update the var
    # names in the solver but really we should not do this
    # TODO: clean up model.sanitize_varnames
    vanilla_model = get_model('optlang-glpk')
    vanilla_model.reactions.EX_glc__D_e.lower_bound = -1 * glc_uptake - glc_uptake_std
    vanilla_model.reactions.EX_glc__D_e.upper_bound = -1 * glc_uptake + glc_uptake_std
    vanilla_model.objective = growth_reaction_id
    fba_sol = vanilla_model.slim_optimize()

    thermo_data, lexicon, compartment_data = get_thermo_data()

    ecoli = ThermoModel(thermo_data=thermo_data,model = vanilla_model,
                              name = name,
                              )

    ecoli.name = name
    ecoli.logger.setLevel(logging.WARNING)
    ecoli.sloppy = True
    # apply_bounds(ecoli,fva)

    ecoli.solver = solver

    annotate_from_lexicon(ecoli, lexicon)
    apply_compartment_data(ecoli, compartment_data)
    # TFA conversion
    ecoli.prepare()
    ecoli.convert(add_displacement = add_displacement)


    ecoli.print_info()
    ecoli.reactions.get_by_id(growth_reaction_id).lower_bound = observed_growth - \
                                                              1*observed_growth_std

    need_relax = False

    ecoli.repair()

    try:
        ecoli.optimize()
    except (AttributeError, SolverError):
        need_relax = True

    if need_relax:
        # final_model, slack_model, relax_table = relax_dgo(ecoli)
        final_model, slack_model, relax_table = relax_dgo(ecoli, in_place = True)
    else:
        final_model = ecoli


    final_model.reactions.get_by_id(growth_reaction_id).lower_bound = 0
    # apply_bounds(ecoli, original_bounds)
    solution = final_model.optimize()
    print('Objective            : {}'.format(final_model.solution.objective_value))
    print(' - Glucose uptake    : {}'.format(final_model.reactions.EX_glc__D_e.flux))

    filepath = 'models/{}'.format(final_model.name)
    save_json_model(final_model, filepath)

    final_model.logger.info('Build complete for model {}'.format(final_model.name))

    return final_model

if __name__ == '__main__':

    tfa_model = create_tfa_model()

    print('Completed')

    from pytfa.analysis.variability import variability_analysis
    tva = variability_analysis(tfa_model)
    tva.to_csv('outputs/tva_iJO1366.csv')

    # Make thermo model
    # make_thermo_model()

    # Save FBA model
    # make_fba_model()
