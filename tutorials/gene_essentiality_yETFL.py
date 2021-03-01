import cobra.io.json
import pytfa.io.json
from etfl.io.json import load_json_model
from etfl.optim.config import standard_solver_config, gene_ko_config
	
import numpy as np
import optlang
	
from cobra.flux_analysis import single_gene_deletion
	
from tqdm import tqdm
	
import time
from gene_essentiality_FBA import complete_medium
	
solver = 'optlang-gurobi'
#
# yeast_fba = cobra.io.json.load_json_model('models/iJO1366_T0E0N0__20180606_121758.json')
# yeast_fba.solver = solver
# yeast_tfa = pytfa.io.json.load_json_model('models/iJO1366_T1E0N0__20180606_121751.json')
# yeast_tfa.solver = solver
	

yeast = load_json_model('models/yeast8_cEFL_2542_enz_128_bins__20200326_152417.json')

complete_medium(yeast)
yeast.solver = solver
	
standard_solver_config(yeast)
yeast.optimize()
max_growth = yeast.solution.objective_value
if max_growth < 0.2: # some problems are difficult to solve
    max_growth = 0.2
    
print('Growth               : {}'.format(yeast.solution.objective_value))
print(' - Ribosomes produced: {}'.format(yeast.solution.raw.EZ_rib))
print(' - RNAP produced: {}'.format(yeast.solution.raw.EZ_rnap))
	
	
## --
def check_feasibility(model):
    '''
    In this function the feasiblity after setting the 10% of max growth for growth
    lower bound is checked.If the problem is infeasible the gene is essential.
    Otherwise, it's non-essential.
    We have different situations:
        1) the feasible solution is found.
        2) the problem is proved to be infeasible.
        3) the time limit reached and a solution has been found.
        4) the time limit reached and a solution has not been found.
    Each case is treated differently.

    Parameters
    ----------
    model : ME-model
        

    Returns
    -------
    growth : nan if problem is infeasible and a found value otherwise.

    '''
    
    try:
        growth = model.slim_optimize()
    except KeyError:
        # I don't know when it happens, but I kept it!
        growth = np.nan
        
    if not np.isnan(growth):
        # the problem is feasible.
        # not to do anything
        return growth
    else:
        if model.solver.problem.Status == 3:
            # the problem is infeasible
            # not to do anything
            return growth
        elif model.solver.problem.Status == 9 and \
            model.solver.problem.SolCount > 0:
            # the time limit reached but a solution has been found.
            growth = 1
            return growth
        elif model.solver.problem.Status == 9 and \
            model.solver.problem.SolCount == 0:
            # the time limit reached and a solution has not been found.
            model.solver.configuration.timeout = 10*3600 # increasing time limit
            growth = model.slim_optimize()
            model.solver.configuration.timeout = 7200 # setting back to the standard value
            return growth
        else:
            # unknown situation
            growth = -1
            return growth
	
def ko_gene(model, gene_id):
    try:
        the_trans = model.get_translation(gene_id)
    except KeyError:
        return None

    
    initial_value = the_trans.upper_bound 
    the_trans.upper_bound = 0
    # We check for the feasibilty of the problem
    growth = check_feasibility(model)
    #
    the_trans.upper_bound = initial_value
    

    return growth


def timeit(method):
    def timed(*args, **kwargs):
        if 'log_data' in kwargs:
            log_data = kwargs.pop('log_data')
        else:
            log_data = dict()

        ti = time.time()
        result = method(*args, **kwargs)
        tf = time.time()


        name = method.__name__
        log_data[name] = tf-ti

        print('{}: {}s'.format(method.__name__, tf-ti))
	
        return result
    return timed
	
@timeit
def ko_etfl(model):
    growth = dict()
    gene_ko_config(model)
    model.growth_reaction.lower_bound = 0.1*max_growth
    model.objective = optlang.symbolics.Zero # Gene essentiality is a feasibility problem


    for g in tqdm(model.genes):
        # with model as this_model:
        ret = ko_gene(model, g.id)
        if ret is None:
            # There is no translation associated to this gene
            model.logger.warning('There is no translation associated '
                                 'to this gene: {}'.format(g.id))
            continue

        growth[g.id] = ret

    model.growth_reaction.lower_bound = 0
    return growth

@timeit
def ko_fba(model):
    return single_gene_deletion(model)

@timeit
def ko_tfa(model):
    return single_gene_deletion(model)


log_data = dict()

# ko_fba(yeast_fba, log_data = log_data)
# ko_tfa(yeast_tfa, log_data = log_data)
values = ko_etfl(yeast, log_data = log_data)

import pandas as pd
df = pd.DataFrame.from_dict(values, orient='index')
df.to_csv('outputs/gene_essentiality_cEFL.csv')