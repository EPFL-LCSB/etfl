from therme.tests.small_model import create_etfl_model
import pandas as pd
import numpy as np
from optlang.exceptions import SolverError
from cobra.util.solver import OptimizationError

from math import log

scales = range(6)
# scales = range(10)

# mu   = pd.DataFrame(index = scales, columns=scales)
# rib  = pd.DataFrame(index = scales, columns=scales)
# rnap = pd.DataFrame(index = scales, columns=scales)
# pgi  = pd.DataFrame(index = scales, columns=scales)

model = dict()
result = dict()
coeffs = dict()
ramp_results = dict()

objs = ['EZ_rib', 'EZ_rnap', 'MR_b4025']
growth = 'Biomass_Ecoli_core'

MAX_MU = 3

def get_active_growth_bounds(model):
    mu = model.solution.x_dict[growth]
    difflist = [abs(mu - x[0]) for x in model.mu_bins]
    min_diff = min(difflist)
    min_ix = difflist.index(min_diff)

    mu_i, (mu_lb, mu_ub) = model.mu_bins[min_ix]

    return mu_i, mu_lb, mu_ub

def glucose_ramp(model, ramp = None):
    if ramp is None:
        ramp = range(0,26)

    glc_rxn = model.reactions.get_by_id('EX_glc__D_e')
    initial_uptake = glc_rxn.lower_bound
    mus=  []
    for this_uptake in ramp:
        glc_rxn.lower_bound = -1*this_uptake

        try:
            sol = model.optimize()
            sol = sol.x_dict[growth]
        except (AttributeError,OptimizationError):
            sol = np.nan

        mus.append(sol)

    return ramp, mus

for o in [growth] + objs:
    result[o] = pd.DataFrame(index = scales, columns=scales)


for k in scales:
    prot_scale = 10**k

    for l in range(k, max(scales)+1):
        mrna_scale = 10**l

        model[k, l] = create_etfl_model(has_thermo=True,
                                        has_neidhardt=False,
                                        mrna_scaling = mrna_scale,
                                        prot_scaling = prot_scale)
        this_model = model[k, l]

        mincoeff = model[k, l].solver.problem.MinCoeff
        maxcoeff = model[k, l].solver.problem.MaxCoeff

        if not hasattr(this_model,'solution'):

            result[growth].loc[k,l] = np.nan
            for o in objs:
                result[o].loc[k,l] = np.nan
                ramp_results[k,l] = glucose_ramp(this_model)
            continue

        coeffs[k,l] = (mincoeff, maxcoeff, np.log10(maxcoeff/mincoeff))

        mu, mu_lb, mu_ub = get_active_growth_bounds(this_model)
        result[growth].loc[k,l] = mu
        # this_model.growth_reaction.lower_bound = mu_lb
        # this_model.growth_reaction.upper_bound = mu_ub
        #
        #
        for o in objs:
            #     this_model.objective = this_model.variables.get(o)
            #     try:
            #         this_model.optimize()
            #         sol = this_model.solution.f
            #     except (AttributeError, SolverError):
            #         sol = np.nan

            if o.startswith('EZ'):
                scaling = prot_scale
            elif o.startswith('MR'):
                scaling = mrna_scale
            else:
                scaling = 1
            val = model[k, l].solution.x_dict[o]
            result[o].loc[k,l] = np.log10( val / scaling) if val >0 else np.nan

        ramp_results[k,l] = glucose_ramp(this_model)


coeffs = pd.DataFrame.from_dict(coeffs, orient='index')
coeffs.columns = ['min', 'max', 'log10(kappa)']

for o,df in result.items():
    df.to_csv('scaling_test/{}.csv'.format(o))

coeffs.to_csv('scaling_test/coeffs.csv')

# Plotting the simulated ranges

import bokeh.plotting as bp
from bokeh.models import Range1d
from bokeh.layouts import gridplot

bp.curdoc().clear()
bp.output_file('scaling_test/ramps.html')
gridsize = len(scales)
grid = []

for k in range(gridsize):
    grid.append([])
    for l in range(gridsize):

        try:
            glc,mu =  ramp_results[k,l]
            p = bp.figure(height = 100, width=100)
            p.line(glc,mu, line_width = 2)

            if k != l:
                p.yaxis.visible = False
                p.xaxis.visible = False

            p.x_range = Range1d(0,glc[-1])
            p.y_range = Range1d(0,MAX_MU)


            grid[k].append(p)

        except KeyError:
            grid[k].append(None)

bp.show(gridplot(grid))




