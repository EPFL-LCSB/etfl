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
from etfl.optim.variables import EnzymeVariable, mRNAVariable
from etfl.optim.constraints import BackwardCatalyticConstraint, \
                                    ForwardCatalyticConstraint
                                    
from pytfa.analysis.chebyshev import chebyshev_center
from etfl.analysis.dynamic import compute_center_generic

from etfl.optim.constraints import ModelConstraint
from pytfa.optim.utils import symbol_sum

from cobra.util.solver import set_objective
from cobra.io.mat import load_matlab_model

cobra_model = load_matlab_model('../organism_data/info_yeast/Y8_3_4_mod_curated.mat')

flux_to_set = 'growth'
# flux_to_set = 'glucose'

solver = 'optlang-gurobi'
#solver = 'optlang-cplex'

GLC_RXN_ID = 'r_1714'
GROWTH_RXN_ID = 'r_4041'

constrainedUptake = ['r_1604','r_1639','r_1873','r_1879','r_1880',
        'r_1881','r_1671','r_1883','r_1757','r_1891','r_1889','r_1810',
        'r_1993','r_1893','r_1897','r_1947','r_1899','r_1900','r_1902',
        'r_1967',
        #'r_1918', Y7 doesn't have a linoleic acid exchange rxn, which
        #was substituted for octadecenoate and octadecynoate for Y5 and Y6,
        #so leave out.
        'r_1903','r_1548','r_1904','r_2028',
        'r_2038','r_1906','r_2067','r_1911','r_1912','r_1913','r_2090',
        'r_1914','r_2106']   # potassium exchange

unconstrainedUptake = ['r_1714','r_1672','r_1654', # ammonium exchange
                    'r_1992', # oxygen exchange
                    'r_2005', # phosphate exchange
                    'r_2060', # sulphate exchange
                    'r_1861', # iron exchange, for test of expanded biomass def
                    'r_1832', # hydrogen exchange
                    'r_2100', # water exchange
                    'r_4593', # chloride exchange
                    'r_4595', # Mn(2+) exchange
                    'r_4596', # Zn(2+) exchange
                    'r_4597', # Mg(2+) exchange
                    'r_2049', # sodium exchange
                    'r_4594', # Cu(2+) exchange
                    'r_4600', # Ca(2+) exchange
                    'r_2020' ]

anyUptake = constrainedUptake + unconstrainedUptake # for objective function

def chemostat_sim(model):      
    # growth media should not be changed, it's minimal mineral
    
    # applying gecko's rxn modifications in order:
    model.reactions.r_1549.upper_bound = 1e-5 # butanediol secretion
    model.reactions.r_2033.upper_bound = 0.05 # pyruvate secretion
    model.reactions.r_1631.upper_bound = 1e-5 # acetaldehyde secretion
    model.reactions.r_1810.upper_bound = 1e-5 # glycine secretion
    model.reactions.r_1634.upper_bound = 0.62 # acetate secretion
    
    # rxn block (gecko)
    model.reactions.r_0659.lower_bound = 0 # isocitrate dehydrogenase (NADP)
    model.reactions.r_0659.upper_bound = 0
    model.reactions.r_2045.lower_bound = 0 # L-serine transport
    
    return

def _va_sim(model):
    model.objective.direction = 'max'
    sol_max = safe_optim(model)

    model.objective.direction = 'min'
    sol_min = safe_optim(model)

    return sol_min, sol_max




def prep_sol(substrate_uptake, model):

    ret = {'obj':model.solution.objective_value,
           'mu':model.solution.fluxes.loc[model.growth_reaction.id],
           'available_substrate':-1*substrate_uptake,
           'uptake':-1*model.solution.fluxes[GLC_RXN_ID]
           }

#    for exch in model.exchanges:
#        ret[exch.id] = model.solution.fluxes.loc[exch.id]
    for rxn in model.reactions:
        ret[rxn.id] = model.solution.fluxes.loc[rxn.id]
    for enz in model.enzymes:
        ret['EZ_'+ enz.id] = model.solution.raw.loc['EZ_'+enz.id]

    return pd.Series(ret)

if __name__ == '__main__':
    # Do things

    
#    uptake_range = pd.Series(np.arange(-15,-16, -1))
    uptake_range = pd.Series(np.arange(-0.5,-15.5,-0.5))
#    growth_range = pd.Series(np.arange(0.025,0.3,0.025)).append(pd.Series(np.arange(0.3,0.4,0.01)))
    growth_range = pd.Series(np.arange(0.025,0.35,0.025)).append(pd.Series(np.arange(0.35,0.42,0.005)))


    model_files = {
#        'cEFL':'yeast8_cEFL_2584_enz_128_bins__20200512_finer.json',
        'vEFL':'yeast8_vEFL_2584_enz_128_bins__20200515_finer.json',
#        'cETFL':'SlackModel yeast8_cETFL_2584_enz_128_bins__20200512_finer.json',
#        'vETFL':'SlackModel yeast8_vETFL_2542_enz_128_bins__20200202_124725.mat.json',   
    }

    models = {k:load_json_model('models/'+v,solver=solver) for k,v in model_files.items()}
    data = {}
    sol = pd.Series()
    chebyshev_variables = None
    BIGM = 1000
    

    for name,model in models.items():
        # growth_uptake_config(model)
        model.warm_start = None
        model.logger.info('Simulating ...')
        start = time()
        if chebyshev_variables is None:
                chebyshev_variables = model.get_variables_of_type(EnzymeVariable) + \
                    model.get_variables_of_type(mRNAVariable)
        chebyshev_center(model, chebyshev_variables,
                     inplace=True,
                     big_m=BIGM,
                     include=[ForwardCatalyticConstraint,BackwardCatalyticConstraint],
                     exclude=[])
        
        if flux_to_set == 'growth':
            tol = 0.01    
            chemostat_sim(model)
            model.reactions.get_by_id(GLC_RXN_ID).upper_bound = 0
            model.reactions.get_by_id(GLC_RXN_ID).lower_bound = -1000
            for gr in growth_range:
                # minimize substrate uptake
                model.objective = symbol_sum([model.reactions.get_by_id(x).reverse_variable \
                                              for x in anyUptake])
                model.objective_direction = 'min'
                
                model.reactions.get_by_id(GROWTH_RXN_ID).lower_bound = gr
                model.reactions.get_by_id(GROWTH_RXN_ID).upper_bound = gr
                temp_sol = safe_optim(model)
                if np.isnan(temp_sol.objective_value):
                    continue
                # fix substrate uptake
                upt = model.objective.value
                expr = model.objective.expression
                sub_cons = model.add_constraint(kind = ModelConstraint, 
                                 hook = model, 
                                 expr = expr,
                                 id_ = 'fix_substrate',
                                 lb = upt - abs(tol * upt),
                                 ub = upt + abs(tol * upt),)
                # fix growth
                fix_growth(model)
                # minimize total sum of fluxes
                model.objective = symbol_sum([model.reactions.get_by_id(x.id).forward_variable + \
                                    model.reactions.get_by_id(x.id).reverse_variable \
                                              for x in cobra_model.reactions \
                                              if x.id != 'r_4050']) # only metabolic reactions and exclude dna reaction as it's not in vETFL
                model.slim_optimize()
                # fix sum of fluxes
                rhs = model.objective.value
                expr = model.objective.expression
                flux_cons = model.add_constraint(kind = ModelConstraint, 
                                 hook = model, 
                                 expr = expr,
                                 id_ = 'fix_tot_flux',
                                 lb = rhs - abs(tol * rhs),
                                 ub = rhs + abs(tol * rhs),)
                # minimize enzyme usage i.e. max dummy enzyme
                obj_expr = symbol_sum([model.enzymes.dummy_enzyme.variable])
                set_objective(model,obj_expr)
                model.objective_direction = 'max'
                
                model.optimize()
                chebyshev_sol = compute_center_generic(model, model.solution) # this also fixes growth with fix_growth function
                new_sol = prep_sol(upt, model)
                sol = pd.concat([sol,new_sol],axis =1)
                # revert the changes
                release_growth(model)
                model.remove_constraint(sub_cons)
                model.remove_constraint(flux_cons)
                
        data[name] = sol
        stop = time()
        print('Elapsed time: {}'.format(stop - start))
        data[name].to_csv('outputs/crabtree_parsi_cheby_{}.csv'.format(name))
