from collections import namedtuple
import pandas as pd
import numpy  as np

from etfl.io.json import load_json_model
from etfl.optim.config import standard_solver_config


from time import time

solver = 'optlang-gurobi'

DefaultSol = namedtuple('DefaultSol', field_names='f')

def get_active_growth_bounds(model):
    mu = model.growth_reaction.flux
    difflist = [abs(mu - x[0]) for x in model.mu_bins]
    min_diff = min(difflist)
    min_ix = difflist.index(min_diff)

    mu_i, (mu_lb, mu_ub) = model.mu_bins[min_ix]

    return mu_i, mu_lb, mu_ub

def safe_optim(model):
    try:
        out = model.optimize()
    except Exception:
        model.logger.warning('Solver status: {}'.format(model.solver.status))
        out = DefaultSol
        out.f = np.nan
    return  out

def _va_sim(model):
    model.objective.direction = 'max'
    sol_max = safe_optim(model)

    model.objective.direction = 'min'
    sol_min = safe_optim(model)

    return sol_min, sol_max


def simulate(available_uptake, model, variables):
    epsilon = model.solver.configuration.tolerances.optimality

    model.logger.info('available_uptake = {}'.format(available_uptake))
    model.reactions.EX_glc__D_e.lower_bound = available_uptake
    model.growth_reaction.lower_bound = 0
    model.growth_reaction.upper_bound = 10

    model.objective = model.growth_reaction.id
    # model.objective = sum(model.objective_variable)
    model.objective.direction = 'max'
    out = safe_optim(model)

    if model.solver.status == 'infeasible':
        ret = {'obj':np.nan,
               'mu': np.nan,
               'mu_lb':np.nan,
               'mu_ub':np.nan,
               'available_substrate':available_uptake,
               'uptake':np.nan,
               'prot_ratio':np.nan,
               'mrna_ratio':np.nan
               }
        for var in variables:
            ret[var + '_lb'] = np.nan
            ret[var + '_ub'] = np.nan
        print('INFEASIBLE SOLUTION AT q={}'.format(available_uptake))
        return pd.Series(ret)

    mu_i, mu_lb, mu_ub = get_active_growth_bounds(model)
    mu = model.growth_reaction.flux

    try:
        prot_ratio = model.interpolation_variable.prot_ggdw.variable.primal
        mrna_ratio = model.interpolation_variable.mrna_ggdw.variable.primal
    except AttributeError:
        # Model without Neidhardt data
        prot_ratio = np.nan
        mrna_ratio = np.nan


    ret = {'obj':model.solution.f,
           'mu': mu,
           'mu_lb':mu_lb,
           'mu_ub':mu_ub,
           'available_substrate':-1*available_uptake,
           'uptake':-1*model.reactions.EX_glc__D_e.flux,
           'prot_ratio':prot_ratio,
           'mrna_ratio':mrna_ratio
           }

    model.growth_reaction.lower_bound = mu - epsilon
    model.growth_reaction.upper_bound = mu + epsilon

    for var in variables:
        model.objective = model.variables.get(var)

        lb, ub = _va_sim(model)

        ret[var + '_lb'] = lb.f
        ret[var + '_ub'] = ub.f

    print(pd.Series(ret))

    return pd.Series(ret)

if __name__ == '__main__':
    # Do things

    variables = [
                'EZ_rib',
                'EZ_rnap',
                # 'EZ_dummy_enzyme',
                # 'MR_dummy_gene',
                 ]


    # uptake_range = pd.Series(np.arange(-1,-40, -1))
    uptake_range = pd.Series(np.arange(-1,-40, -1))

    models = {
        'T1E1N0': load_json_model('models/RelaxedModel iJO1366_T1E1N0_431_enz_128_bins__20180926_140929.json',
                                  solver=solver),
        'T1E1N1': load_json_model('models/RelaxedModel iJO1366_T1E1N1_431_enz_128_bins__20180926_124941.json',
                                  solver=solver),
        'T0E1N0':load_json_model('models/iJO1366_T0E1N0_431_enz_128_bins__20180926_144307.json',
                                 solver = solver),
        'T0E1N1': load_json_model('models/iJO1366_T0E1N1_431_enz_128_bins__20180926_135704.json',
                                   solver = solver),

       # 'T0E1N0':load_json_model('models/iJO1366_T0E1N0_430_enz_256_bins__20180903_130235.json',
       #                          solver = solver),
       # 'T1E1N0':load_json_model('models/RelaxedModel '
       #                          'iJO1366_T1E1N0_430_enz_256_bins__20180903_131226.json',
       #                          solver = solver),
       #  # 'T0E1N1':load_json_model('models/iJO1366_T0E1N1_431_enz_256_bins__20180903_152226.json',
       #  #                           solver = solver),
       # 'T1E1N1':load_json_model('models/RelaxedModel '
       #                          'iJO1366_T1E1N1_430_enz_256_bins__20180903_141648.json',
       #                          solver = solver),
        # 'T0E1N0':load_json_model('models/iJO1366_T0E1N0_350_enz_256_bins__20180830_112731.json',
        #                          solver = solver),
        # 'T1E1N0':load_json_model('models/RelaxedModel '
        #                          'iJO1366_T1E1N0_350_enz_256_bins__20180830_113554.json',
        #                          solver = solver),
        # 'T0E1N1':load_json_model('models/iJO1366_T0E1N1_350_enz_256_bins__20180830_111752.json',
        #                          solver = solver),
        # 'T1E1N1':load_json_model('models/RelaxedModel '
        #                          'iJO1366_T1E1N1_350_enz_256_bins__20180830_121200.json',
        #                          solver = solver),
              }
    data = {}
    for name,model in models.items():

        standard_solver_config(model)

        model.logger.info('Simulating ...')
        start = time()
        data[name] = uptake_range.apply(simulate, args=[model,variables])
        stop = time()
        print('Elapsed time: {}'.format(stop - start))
        data[name].to_csv('outputs/benchmark_{}.csv'.format(name))
