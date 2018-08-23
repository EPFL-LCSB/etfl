from therme.tests.small_model import create_etfl_model
from therme.analysis.dynamic import run_dynamic_etfl
from therme.optim.config import standard_solver_config

import bokeh.plotting as bp
from bokeh.palettes import Category10
from bokeh.layouts import row, gridplot

# model = create_etfl_model(  1,0,
#                             prot_scaling=1e3,
#                             mrna_scaling=1e3,
#                             n_mu_bins = 128,
#                             )
# model.reactions.ATPM.lower_bound = 4
# model.reactions.EX_glc__D_e.lower_bound = -1
# model.optimize()
# model.objective = model.growth_reaction.forward_variable + \
#                   model.reactions.EX_glc__D_e.reverse_variable


from therme.io.json import load_json_model

model = load_json_model('models/iJO1366_T0E1N0_473_enz_256_bins__20180822_081415.json')
# model = load_json_model('models/iJO1366_T0E1N0_350_enz_256_bins__20180802_123525.json')
# model = load_json_model('models/iJO1366_T0E1N1_350_enz_256_bins__20180802_124700.json')

def plot_dynamics(model, time_data):

    bp.curdoc().clear()

    p1 = bp.figure()
    p2 = bp.figure()
    p3 = bp.figure(y_axis_type = 'log')
    p4 = bp.figure(y_axis_type = 'log')

    t = time_data.loc['t']

    vars_to_obs = ['MR_b4025', 'MR_b1779', 'EZ_PGI_CPLX0_7877', 'EZ_GAPD_GAPDH_A_CPLX', 'EZ_rib', 'EZ_rnap']

    cmap = Category10[10]

    for e, the_var in enumerate(['S_EX_glc__D_e','S_EX_glyc_e','X']):
        plot_hist(p1, t, time_data, the_var, cmap[e])

    rev_glc = model.reactions.EX_glc__D_e.reverse_variable.name
    for e, the_var in enumerate(['PGI','GAPD',
                                 'ub_EX_glc__D_e', 'ub_EX_glyc_e',
                                 'act_EX_glc__D_e', 'act_EX_glyc_e',
                                 'mu']):
        plot_hist(p2, t, time_data, the_var, cmap[e])
    for e, the_var in enumerate(vars_to_obs):
        plot_hist(p3, t, time_data, the_var, cmap[e])

    effluxes = time_data.loc[ time_data.index.str.startswith('EX_')
                            & time_data.index.str.endswith('_e') ]

    effluxes['score'] = effluxes.std(axis=1) / effluxes.mean(axis=1)

    effluxes = effluxes.reindex(
        effluxes['score'].sort_values(ascending = False).index,
        axis=0)

    for e,(the_var, values) in enumerate(effluxes[:5].T.iteritems()):

        p4.line(t, values,
                color = cmap[e],
                line_width = 2,
                legend=the_var)

    for p in [p1,p2,p3,p4]:
        p.legend.location = 'top_left'

    bp.output_file('plots/dETFL_small.html')
    bp.show(gridplot([[p1,p2],[p3,p4]]))


def plot_hist(figure, t, time_data, the_var, color):
    figure.line(t, time_data.loc[the_var],
                color=color,
                line_width=2,
                legend=the_var)
    figure.circle(t, time_data.loc[the_var],
                  color=color,
                  legend=the_var)


def get_uptake_funs():

    uptake_funs = dict()

    # Glucose:

    Vmax0 = 10 # mmol/(h.mmol[E]) Mahadevan et al. 2002
    # Vmax = 5 # For testing
    Km0 = 0.015 # mM, Mahadevan et al. 2002, Wong et al. 1997
    # Km = 1.5 # For testing

    uptake_funs['EX_glc__D_e'] = lambda x: Vmax0 * x / (Km0 + x)

    # Glycerol:

    Vmax1 = 10*60 / 1000 * 56106 # umol/(min.mg[E]) --> mmol/(h.mmol[E])
                                # Voegele et al. 1992
                                # Pettigrew et al. 1987
    Vmax1 = 10
    Km1 = 0.010 # mM, Voegele et al. 1992

    uptake_funs['EX_glyc_e'] = lambda x: Vmax1 * x / (Km1 + x)

    return uptake_funs



if __name__ == '__main__':

    timestep = 0.1

    S0 = 2.0 #mmol/L
    S1 = 10.0 #mmol/L
    # X0 = 0.1 #mmol[Cell]/L
    X0 = 0.001 #g[Cell]/L

    glc_step = lambda t,S: S+S1 if abs(t-1)<=timestep and S<=S0 else max(S,0)
    gly_step = lambda t,S: S0 if (t-1)<=timestep and S<=S0 else max(S,0)

    glyc_glc_switch = {
        'EX_glc__D_e': glc_step,
        'EX_glyc_e'  : gly_step
    }

    standard_solver_config(model)

    time_data = run_dynamic_etfl(model,
                                 timestep=timestep,
                                 tfinal=5,
                                 uptake_fun=get_uptake_funs(),
                                 medium_fun = glyc_glc_switch,
                                 S0={'EX_glyc_e':S0, 'EX_glc__D_e':0.0},
                                 X0=X0,
                                 inplace=True,
                                 )

    plot_dynamics(model, time_data)
