from etfl.tests.small_model import create_etfl_model
from etfl.optim.config import standard_solver_config
from pytfa.optim.debugging import find_maxed_vars

m = create_etfl_model(0,1)

mv = find_maxed_vars(m)
from etfl.io.dict import model_to_dict, model_from_dict

dm = model_to_dict(m)

md = model_from_dict(dm)

standard_solver_config(md)

md.optimize()

print(abs(m.solution.f - md.solution.f))

print('Objective            : {}'.format(md.solution.f))
print(' - Glucose uptake    : {}'.format(md.reactions.EX_glc__D_e.flux))
print(' - Growth            : {}'.format(md.growth_reaction.flux))
print(' - Ribosomes produced: {}'.format(md.ribosome.X))
print(' - RNAP produced: {}'.format(md.rnap.X))


print(mv)


assert abs(m.solution.f - md.solution.f) < md.solver.configuration.tolerances.optimality