from etfl.tests.small_model import create_simple_dynamic_model
from etfl.analysis.dynamic import run_dynamic_etfl

smol_model = create_simple_dynamic_model()

def set_mix(glc_amount, lcts_amount):
    smol_model.reactions.EX_glc.bounds  = (-glc_amount ,0)
    smol_model.reactions.EX_lcts.bounds = (-lcts_amount,0)

def get_yield(sol):
    return sol.objective_value / (sol.fluxes.glc2mint / 6 + sol.fluxes.lcts2mint / 12)

set_mix(10,0)
sol_g = smol_model.optimize()

set_mix(0,5)
sol_l = smol_model.optimize()

# Yield with glc, yield with lcts
Yg = get_yield(sol_g)
Yl = get_yield(sol_l)

# Initial solution on lcts/glc mix
set_mix(10,5)
ini_sol = smol_model.optimize()

uptake_fun = {'EX_glc' :lambda x:10, # No point in putting kinetics, we just put fixed LBs
              'EX_lcts':lambda x:5}
medium_fun = {'EX_glc' :lambda t,x:max(0,x), # no negative concentrations
              'EX_lcts':lambda t,x:max(0,x)}

S0 = {'EX_glc' : 1, # No point in putting kinetics, we just put fixed LBs
      'EX_lcts': 1}
X0 = 1

# Dyno:
tsol =  run_dynamic_etfl(smol_model, timestep=0.01, tfinal=1, uptake_fun=uptake_fun, medium_fun=medium_fun,
                         uptake_enz=[],
                         S0=S0, X0=X0,
                         inplace=True,
                         initial_solution = ini_sol,
                         chebyshev_include=False,
                         dynamic_constraints={
                                'mRNA_degradation':False,
                                'mRNA_synthesis':False,
                                'enzyme_synthesis':True,
                                'enzyme_degradation':True,
                            },
                         mode='forward')

tsol.to_csv('outputs/tsol_small_model.csv')
