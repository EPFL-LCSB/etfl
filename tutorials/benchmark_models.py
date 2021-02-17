from collections import namedtuple
import pandas as pd
import numpy  as np

from etfl.io.json import load_json_model
from etfl.optim.config import standard_solver_config, growth_uptake_config

from etfl.optim.variables import GrowthActivation, BinaryActivator

from time import time
from copy import copy

from etfl.optim.utils import fix_growth, release_growth, \
                            get_active_growth_bounds, safe_optim

try:
    from gurobipy import GRB
except ModuleNotFoundError:
    pass

solver = 'optlang-gurobi'
# solver = 'optlang-cplex'


def _va_sim(model):
    model.objective.direction = 'max'
    sol_max = safe_optim(model)

    model.objective.direction = 'min'
    sol_min = safe_optim(model)

    return sol_min, sol_max


def simulate(available_uptake, model, variables, warm_start=None):

    model.logger.info('available_uptake = {}'.format(available_uptake))
    model.reactions.EX_glc__D_e.lower_bound = available_uptake
    model.growth_reaction.lower_bound = 0
    model.growth_reaction.upper_bound = 10

    model.objective = model.growth_reaction.id
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

    growth_solution = copy(model.solution)
    mu_i, mu_lb, mu_ub = get_active_growth_bounds(model)
    mu = model.growth_reaction.flux
    # release_warm_start(model)

    try:
        prot_ratio = model.interpolation_variable.prot_ggdw.variable.primal
        mrna_ratio = model.interpolation_variable.mrna_ggdw.variable.primal
    except AttributeError:
        # Model without Neidhardt data
        prot_ratio = np.nan
        mrna_ratio = np.nan


    ret = {'obj':model.solution.objective_value,
           'mu': mu,
           'mu_lb':mu_lb,
           'mu_ub':mu_ub,
           'available_substrate':-1*available_uptake,
           'uptake':-1*growth_solution.fluxes['EX_glc__D_e'],
           'prot_ratio':prot_ratio,
           'mrna_ratio':mrna_ratio
           }

    fix_growth(model, model.solution)

    for var in variables:
        model.objective = model.variables.get(var)

        lb, ub = _va_sim(model)

        ret[var + '_lb'] = lb.objective_value
        ret[var + '_ub'] = ub.objective_value

    print(pd.Series(ret))

    release_growth(model)
    # apply_warm_start(model, growth_solution)

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
    # uptake_range = pd.Series(np.arange(-1,-30, -1))
    uptake_range = pd.Series(np.arange(-1,-20, -1))

    model_files = {
        # 'EFL':'iJO1366_EFL_1783_enz_128_bins__20190221_165945.json',
        # 'ETFL':'SlackModel iJO1366_ETFL_1783_enz_128_bins__20190221_172509.json',
        # 'vEFL':'iJO1366_vEFL_1783_enz_128_bins__20190221_182321.json',
        # 'vETFL':'SlackModel iJO1366_vETFL_1783_enz_128_bins__20190221_185154.json',
        # 'vETFL65':'SlackModel iJO1366_vETFL_1783_enz_128_bins__20190221_194424.json',
        # 'vETFL_infer':'SlackModel iJO1366_vETFL_infer_2088_enz_128_bins__20190221_204202.json',
        # 'vETFL65_infer': 'SlackModel iJO1366_vETFL_infer_2088_enz_128_bins__20190221_220847.json',
        #
        'EFL':'iJO1366_EFL_431_enz_128_bins__20190121_080047.json',
        'ETFL':'SlackModel iJO1366_ETFL_431_enz_128_bins__20190121_110027.json',
        'vEFL':'iJO1366_vEFL_431_enz_128_bins__20190121_090316.json',
        'vETFL':'SlackModel iJO1366_vETFL_431_enz_128_bins__20190122_145700.json',
        'vETFL65':'SlackModel iJO1366_vETFL_431_enz_128_bins__20190122_155755.json',
        'vETFL_infer':'SlackModel iJO1366_vETFL_2084_enz_128_bins__20190122_170118.json',
        'vETFL65_infer':'SlackModel iJO1366_vETFL_2084_enz_128_bins__20190124_082824.json',
    }

    models = {k:load_json_model('models/'+v,solver=solver) for k,v in model_files.items()}
    data = {}

    for name,model in models.items():
        growth_uptake_config(model)
        model.warm_start = None
        model.logger.info('Simulating ...')
        start = time()
        data[name] = uptake_range.apply(simulate, args=[model,variables])
        stop = time()
        print('Elapsed time: {}'.format(stop - start))
        data[name].to_csv('outputs/benchmark_{}.csv'.format(name))
