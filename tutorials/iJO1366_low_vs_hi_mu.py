from therme.io.json import load_json_model

from therme.optim.variables import mRNAVariable, EnzymeVariable
from therme.optim.utils import fix_integers

from pytfa.analysis import  variability_analysis,           \
                            apply_reaction_variability,     \
                            apply_generic_variability,      \
                            apply_directionality

from therme.analysis.utils import enzymes_to_peptides_conc

import pandas as pd


ecoli = load_json_model('models/RelaxedModel iJO1366_T1E1N1_431_enz_128_bins__20180926_124941.json')
ecoli.optimize()

def print_sol(model):
    print('Objective            : {}'.format(model.solution.f))
    print(' - Glucose uptake    : {}'.format(model.solution.x_dict.EX_glc__D_e))
    print(' - Growth            : {}'.format(model.solution.x_dict[model.growth_reaction.id]))
    print(' - Ribosomes produced: {}'.format(model.solution.raw.EZ_rib))
    print(' - RNAP produced: {}'.format(model.solution.raw.EZ_rnap))
    try:
        print(' - DNA produced: {}'.format(model.solution.raw.DN_DNA))
    except AttributeError:
        pass

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

    print_sol(continuous_model)

    variables = EnzymeVariable

    eva = variability_analysis(continuous_model, variables)
    peptides_conc_min = pd.Series(enzymes_to_peptides_conc(continuous_model, eva['minimum']))
    peptides_conc_max = pd.Series(enzymes_to_peptides_conc(continuous_model, eva['maximum']))
    peptides_conc = pd.concat([peptides_conc_min,peptides_conc_max], axis=1)
    peptides_conc.to_csv('outputs/iJO_T1E1N1_low_hi_{}_pep.csv'.format(glc_uptake))

    rescale = lambda row: ecoli.mrnas.get_by_id(row.name[3:]).scaling_factor * row
    eva_real = eva.apply(rescale, axis=1)
    eva_real.to_csv('outputs/iJO_T1E1N1_low_hi_{}_enz.csv'.format(glc_uptake))

    # mva = variability_analysis(continuous_model, mRNAVariable)
    # rescale = lambda row: ecoli.mrnas.get_by_id(row.name[3:]).scaling_factor * row
    # mva_real = mva.apply(rescale, axis=1)
    # mva_real.to_csv('outputs/iJO_T1E1N1_low_hi_{}_mrna.csv'.format(glc_uptake))




