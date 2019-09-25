import cobra
from etfl.io.json import load_json_model
from etfl.optim.config import standard_solver_config

from etfl.core.reactions import TranscriptionReaction, TranslationReaction, \
    ProteinComplexation, DegradationReaction

from cobra.flux_analysis import flux_variability_analysis
from etfl.optim.utils import fix_growth, release_growth, \
                            get_active_growth_bounds, safe_optim

from tqdm import tqdm

from copy import copy
from numpy import isnan
import pandas as pd

solver = 'optlang-gurobi'

EPSILON = 1e-6

# vETFL
ecoli = load_json_model('models/SlackModel iJO1366_vETFL_431_enz_128_bins__20190701_082518.json',
                        solver = solver)

expression_reaction_classes = [ProteinComplexation, DegradationReaction,
                               TranscriptionReaction, TranslationReaction]

iJO1366 = cobra.io.load_json_model('iJO1366_with_xrefs.json')
fva_fba = flux_variability_analysis(iJO1366)
reaction_list = fva_fba[fva_fba['minimum'] * fva_fba['maximum'] < 0].index


# reaction_list = [r for r in ecoli.reactions
#                  if r.lower_bound*r.upper_bound < 0 and \
#                  not any([isinstance(r,cls)
#                              for cls in expression_reaction_classes])]

print('* VA will be performed on',len(reaction_list),'reactions.')


standard_solver_config(ecoli, verbose=False)
ecoli.optimize()

print('Objective            : {}'.format(ecoli.solution.objective_value))
print(' - Glucose uptake    : {}'.format(ecoli.reactions.EX_glc__D_e.flux))
print(' - Growth            : {}'.format(ecoli.growth_reaction.flux))
print(' - Ribosomes produced: {}'.format(ecoli.ribosome.X))
print(' - RNAP produced: {}'.format(ecoli.rnap.X))

growth_solution = copy(ecoli.solution)
mu_i, mu_lb, mu_ub = get_active_growth_bounds(ecoli)
mu = ecoli.growth_reaction.flux

def _va_sim(model, rxn, epsilon = EPSILON):
    lb, ub = rxn.bounds

    rxn.lower_bound = epsilon
    # model.objective.direction = 'max'
    sol_max = safe_optim(model)
    rxn.lower_bound = lb

    rxn.upper_bound = -epsilon
    # model.objective.direction = 'min'
    sol_min = safe_optim(model)
    rxn.upper_bound = ub

    return ~isnan(sol_min.objective_value), \
           ~isnan(sol_max.objective_value)

fix_growth(ecoli, ecoli.solution)

results = dict()
ecoli.objective = 0
for rxnid in tqdm(reaction_list):
    rxn = ecoli.reactions.get_by_id(rxnid)#.id)
    lb, ub = _va_sim(ecoli, rxn)
    results[rxn.id] = (lb,ub)

release_growth(ecoli)

#fva = flux_variability_analysis(model=ecoli, reaction_list=reaction_list)
fva = pd.DataFrame.from_dict(results, orient = 'index')
fva.to_csv('outputs/fva_etfl.csv')