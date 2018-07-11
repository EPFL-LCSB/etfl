from therme.io.json import load_json_model

from therme.optim.variables import mRNAVariable
from therme.optim.utils import fix_integers
from pytfa.analysis import  variability_analysis,           \
                            apply_reaction_variability,     \
                            apply_generic_variability,       \
                            apply_directionality



ecoli = load_json_model('models/RelaxedModel iJO1366_T1E1N1_346_enz_256_bins__20180710_071309.json')
ecoli.optimize()

print('Growth               : {}'.format(ecoli.solution.f))
print(' - Ribosomes produced: {}'.format(ecoli.solution.x_dict.EZ_rib))
print(' - RNAP produced: {}'.format(ecoli.solution.x_dict.EZ_rnap))

# ID	    109686
# Property	Glucose uptake rate of strain C-3000 in minimal M9 media
# Organism	Bacteria Escherichia coli
# Value	    12
# Range	    ±0.5
# Units	    mmol/gDW/h
# Reference	Jain R, Srivastava R. Metabolic investigation of host/pathogen interaction using MS2-infected Escherichia coli. BMC Syst Biol. 2009 Dec 30 3 :121. doi: 10.1186/1752-0509-3-121. p.5 right column 3rd paragraph
# Reference PubMed ID	20042079
# Method	"Changes in the glucose concentration of the medium were used to determine the glucose uptake rate of uninfected Escherichia coli C-3000 cells in minimal M9 media."
# Comments	"The glucose uptake rate of the uninfected cells was 12±0.5mmol/(gDW?h) (n=3) and the oxygen uptake rate was 27±1mmol/(gDW?h) (n=2)."
uptake_high = -12.5 #mmol/gDW/h

uptake_low  = -1 # mmol/gDW/h

# ecoli.growth_reaction.lower_bound = ecoli.solution.f - ecoli.solver.configuration.tolerances.feasibility
ecoli.optimize()


ecoli.solver.configuration.verbosity = 1
ecoli.solver.configuration.tolerances.feasibility = 1e-9
try:
    ecoli.solver.problem.Params.NumericFocus = 3
except AttributeError:
    pass

ecoli.solver.configuration.presolve = True

# print('Fixing integers..')
# continuous_model = fix_integers(ecoli)
continuous_model = ecoli

# continuous_model.remove_constraint(continuous_model.sos1_constraint.interpolation_integer_SOS1)
print('optimizing..')
continuous_model.optimize()


def get_mu_bin(model, mu):
    closeness = [abs(x[0]-mu) for x in ecoli.mu_bins]
    ix = closeness.index(min(closeness))
    return model.mu_bins[ix][1]

# Calculate variability analysis on all continuous variables

for glc_uptake in [uptake_high,uptake_low]:
    continuous_model.growth_reaction.lower_bound = 0
    continuous_model.growth_reaction.upper_bound = 10
    continuous_model.reactions.EX_glc__D_e.lower_bound = glc_uptake - 0.1
    continuous_model.reactions.EX_glc__D_e.upper_bound = glc_uptake + 0.1

    continuous_model.optimize()
    mu = continuous_model.growth_reaction.flux
    mu_lo, mu_hi = get_mu_bin(continuous_model, mu)

    continuous_model.growth_reaction.upper_bound = mu_hi
    continuous_model.growth_reaction.lower_bound = mu_lo

    print('Growth               : {}'.format(continuous_model.solution.f))
    print(' - Ribosomes produced: {}'.format(
        continuous_model.solution.x_dict.EZ_rib))
    print(
        ' - RNAP produced: {}'.format(continuous_model.solution.x_dict.EZ_rnap))

    variables = mRNAVariable

    eva = variability_analysis(continuous_model, variables)
    eva.to_csv('outputs/iJO_T1E1N1_low_hi_{}.csv'.format(glc_uptake))

    # rva = variability_analysis(continuous_model, 'reactions')
    # rva.to_csv('outputs/iJO_T1E1N1_low_hi_{}_rxns.csv'.format(glc_uptake))




