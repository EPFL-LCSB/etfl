from therme.tests.small_model import create_etfl_model
from therme.analysis.dynamic import run_dynamic_etfl

model = create_etfl_model(  0,0,
                            prot_scaling=1e0,
                            mrna_scaling=1e2,)

model.reactions.EX_glc__D_e.lower_bound = -1
model.optimize()

time_data = run_dynamic_etfl(model,
                             initial_solution=model.solution.x_dict,
                             timestep = 0.1,
                             tfinal = 1
                            )


