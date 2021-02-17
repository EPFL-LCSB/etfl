from collections import namedtuple
import pandas as pd
import numpy  as np

from etfl.io.json import load_json_model
from etfl.optim.config import standard_solver_config, growth_uptake_config

from etfl.optim.variables import GrowthActivation, BinaryActivator

from pytfa.optim.utils import symbol_sum


from time import time
from copy import copy

from etfl.optim.utils import fix_growth, release_growth, \
                            get_active_growth_bounds, safe_optim
                            
from etfl.optim.variables import ModelVariable,EnzymeVariable, mRNAVariable
from etfl.optim.constraints import ModelConstraint

try:
    from gurobipy import GRB
except ModuleNotFoundError:
    pass

solver = 'optlang-gurobi'
# solver = 'optlang-cplex'

class TotalResourse(ModelVariable):
    """
    Represents a variable for total RNA or protein
    """
    
    prefix = 'TOT_'
    
class TotalResourseConstraint(ModelConstraint):
    """
    Represents a variable for total RNA or protein
    """
    
    prefix = 'TTC_'

def _va_sim(model):
    model.objective.direction = 'max'
    sol_max = safe_optim(model)

    model.objective.direction = 'min'
    sol_min = safe_optim(model)

    return sol_min, sol_max


def simulate(available_uptake, model, variables, warm_start=None):

#    model.solver.problem.reset()
    model.logger.info('available_uptake = {}'.format(available_uptake))
    model.reactions.r_1714.lower_bound = available_uptake
    model.reactions.r_1714.upper_bound = available_uptake
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
#    mu_i, mu_lb, mu_ub = get_active_growth_bounds(model)
    mu = model.growth_reaction.flux
    # release_warm_start(model)

    try:
        prot_ratio = model.interpolation_variable.prot_ggdw.variable.primal
        mrna_ratio = model.interpolation_variable.mrna_ggdw.variable.primal
        dna_ratio = model.interpolation_variable.dna_ggdw.variable.primal
        lipid_ratio = model.interpolation_variable.lipid_ggdw.variable.primal
        carbohydrate_ratio = model.interpolation_variable.carbohydrate_ggdw.variable.primal
        ion_ratio = model.interpolation_variable.ion_ggdw.variable.primal
    except AttributeError:
        # Model without Neidhardt data
        prot_ratio = model.variables.TOT_prot.primal
        mrna_ratio = model.variables.TOT_RNA.primal
        dna_ratio = np.nan
        lipid_ratio = np.nan
        carbohydrate_ratio = np.nan
        ion_ratio = np.nan

    ret = {'obj':model.solution.objective_value,
           'mu': mu,
#           'mu_lb':mu_lb,
#           'mu_ub':mu_ub,
           'available_substrate':-1*available_uptake,
           'uptake':-1*growth_solution.fluxes['r_1714'],
           'prot_ratio':prot_ratio,
           'mrna_ratio':mrna_ratio,
           'dna_ratio':dna_ratio,
           'carbohydrate_ratio':carbohydrate_ratio,
           'lipid_ratio':lipid_ratio,
           'ion_ratio':ion_ratio,
           }

    fix_growth(model, model.solution)

    for var in variables:
        # THIS WILL DO THE VA ON ETHANOL
        model.objective = model.variables.get(var)

        lb, ub = _va_sim(model)

        ret[var + '_lb'] = lb.objective_value#np.nan
        ret[var + '_ub'] = ub.objective_value#np.nan

    print(pd.Series(ret))

    release_growth(model)
    # apply_warm_start(model, growth_solution)
    
    # Add values of other secretions in the ret dictionnary
    for rxn in model.reactions:
        ret[rxn.id] = model.solution.fluxes.loc[rxn.id]
    for enz in model.enzymes:
        ret['EZ_'+ enz.id] = growth_solution.raw.loc['EZ_'+ enz.id]
    for mRNA in model.mrnas:
        ret['MR_'+ mRNA.id] = growth_solution.raw.loc['MR_'+ mRNA.id]

    return pd.Series(ret)

def VA_prepare(model):
    # prepare the model for mRNA and protein VA
    Enz_vars = model.get_variables_of_type(EnzymeVariable)
    total_prot = model.add_variable(kind = TotalResourse,
                                   hook = model,
                                   id_ = 'prot',
                                   lb = 0,
                                   ub = 1)
    
    expr = symbol_sum([x for x in Enz_vars])
    model.add_constraint(kind = TotalResourseConstraint, 
                                 hook = model, 
                                 expr = expr - total_prot,
                                 id_ = 'prot',
                                 lb = 0,
                                 ub = 0)
    
    RNA_vars = model.get_variables_of_type(mRNAVariable)
    total_rna = model.add_variable(kind = TotalResourse,
                                   hook = model,
                                   id_ = 'RNA',
                                   lb = 0,
                                   ub = 1)
    
    expr = symbol_sum([x for x in RNA_vars])
    model.add_constraint(kind = TotalResourseConstraint, 
                                 hook = model, 
                                 expr = expr - total_rna,
                                 id_ = 'RNA',
                                 lb = 0,
                                 ub = 0)
    

if __name__ == '__main__':
    # Do things

    variables = [
#                'EZ_rib',
#                'EZ_rnap',
#                'r_1761',    # ethanol exchange
#                'MR_dummy_gene',
#                'EZ_dummy_enzyme',
                 ]


    # uptake_range = pd.Series(np.arange(-15,-20, -1))
    uptake_range = pd.Series(np.arange(-1/3,-5,-1/3)).append(pd.Series(np.arange(-5,-16,-1)))
    # uptake_range = pd.Series(np.arange(-10,-11,-1))

    model_files = {
#        'cEFL':'yeast8_cEFL_2542_enz_128_bins__20200326_152417.json',
        'cETFL':'SlackModel yeast8_cETFL_2542_enz_128_bins__20200326_152515.json',
#        'vEFL':'yeast8_vEFL_2542_enz_128_bins__20200326_191004.json',
        'vETFL':'SlackModel yeast8_vETFL_2542_enz_128_bins__20200326_190837.json',
    }

    models = {k:load_json_model('models/'+v,solver=solver) for k,v in model_files.items()}
    data = {}

    for name,model in models.items():
        # growth_uptake_config(model)
        model.warm_start = None
        model.logger.info('Simulating ...')
        start = time()
        # Add a new variable for total mRNA and total protein to do VA
        VA_prepare(model)
        data[name] = uptake_range.apply(simulate, args=[model,variables])
        stop = time()
        print('Elapsed time: {}'.format(stop - start))
        data[name].to_csv('outputs/benchmark_{}.csv'.format(name))
