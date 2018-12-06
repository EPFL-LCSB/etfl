from etfl.tests.small_model import create_etfl_model
from etfl.analysis.dynamic import run_dynamic_etfl
from etfl.optim.config import standard_solver_config
from etfl.core.reactions import EnzymaticReaction
from etfl.optim.constraints import ForwardCatalyticConstraint, \
    BackwardCatalyticConstraint, ExpressionCoupling

import bokeh.plotting as bp
from bokeh.palettes import Category10
from bokeh.layouts import row, gridplot

from bokeh.io import export_svgs

import pandas as pd
#
# model = create_etfl_model(  0,0,
#                             n_mu_bins = 128,
#                             )
# # model.reactions.ATPM.lower_bound = 4
# model.reactions.EX_glc__D_e.lower_bound = -1
# model.optimize()
# model.objective = model.growth_reaction.forward_variable + \
#                   model.reactions.EX_glc__D_e.reverse_variable


from etfl.io.json import load_json_model

model = load_json_model('models/iJO1366_T0E1N1_431_enz_128_bins__20180926_135704.json')
# model = load_json_model('models/RelaxedModel iJO1366_T1E1N1_431_enz_128_bins__20180926_124941.json')

def get_enzymes_of_subsystem(model, subsystem):
    reactions = [x for x in model.reactions if subsystem.lower() in x.subsystem.lower()]

    enzymes = [item for rxn in reactions if hasattr(rxn,'enzymes')
               for item in rxn.enzymes
               ]

    return enzymes

def get_mrna_total(time_data, mrnas):
    return time_data.loc[[x.variable.name * x.scaling_factor
                          for x in mrnas]].sum()

def get_prot_total(time_data, enzymes):
    enz = time_data.loc[[x.variable.name for x in enzymes]]
    sigma =  pd.Series([x.scaling_factor for x in enzymes], index = enz.index)

    return (sigma*enz.T).T.sum()


def plot_dynamics(model, time_data):

    bp.curdoc().clear()

    # time_data += model.solver.configuration.tolerances.feasibility

    p1 = bp.figure(y_axis_type = 'log')
    p2 = bp.figure(y_axis_type = 'log')
    p3 = bp.figure(y_axis_type = 'log')
    p4 = bp.figure(y_axis_type = 'log')
    pex = bp.figure()

    t = time_data.loc['t']

    vars_to_obs = ['MR_b4025', 'MR_b1779', 'EZ_rib', 'EZ_rnap']

    cmap = Category10[10]

    for e, the_var in enumerate(['S_EX_glc__D_e','S_EX_glyc_e','X']):
        plot_hist(p1, t, time_data, the_var, cmap[e])

    rev_glc = model.reactions.EX_glc__D_e.reverse_variable.name
    for e, the_var in enumerate(['ub_EX_glc__D_e', 'ub_EX_glyc_e',
                                 'act_EX_glc__D_e', 'act_EX_glyc_e',
                                 'mu']):
        plot_hist(p2, t, time_data, the_var, cmap[e])

    glycolysis_enzymes = get_enzymes_of_subsystem(model, 'Glycolysis')
    ppp_enzymes = get_enzymes_of_subsystem(model, 'pentose phosphate')
    tca_enzymes = get_enzymes_of_subsystem(model, 'Citric Acid Cycle')

    for e, the_var in enumerate(['ub_EX_glc__D_e', 'ub_EX_glyc_e',
                                 'act_EX_glc__D_e', 'act_EX_glyc_e',
                                 'mu']):
        plot_hist(p2, t, time_data, the_var, cmap[e])

    for e, the_var in enumerate(vars_to_obs):
        plot_hist(p3, t, time_data, the_var, cmap[e])

    time_data.loc['glycolysis'] = get_prot_total(time_data, glycolysis_enzymes)
    time_data.loc['ppp'] = get_prot_total(time_data, ppp_enzymes)
    time_data.loc['tca'] = get_prot_total(time_data, tca_enzymes)
    plot_hist(p4, t, time_data, 'glycolysis', cmap[e+1])
    plot_hist(p4, t, time_data, 'ppp', cmap[e+2])
    plot_hist(p4, t, time_data, 'tca', cmap[e+3])

    effluxes = time_data.loc[ time_data.index.str.startswith('EX_')
                              & time_data.index.str.endswith('_e') ]

    effluxes['score'] = effluxes.std(axis=1) / effluxes.mean(axis=1)

    effluxes = effluxes.reindex(
        effluxes['score'].sort_values(ascending = False).index,
        axis=0)

    for e,(the_var, values) in enumerate(effluxes[:5].T.iteritems()):

        pex.line(t, values,
                 color = cmap[e],
                 line_width = 2,
                 legend=the_var)

    for p in [p1,p2,p3,p4,pex]:
        #     p.legend.location = 'top_left'
        p.output_backend = 'svg'

    gp = gridplot([[p1,p3],[p2,p4],[pex]])
    bp.output_file('plots/dETFL_cheby.html')
    bp.show(gp)
    export_svgs(gp, filename='plots/dETFL_cheby.svg')


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
    # Vmax0 = 15 # mmol/(h.mmol[E]) for testing
    # Vmax = 5 # For testing
    Km0 = 0.015 # mM, Mahadevan et al. 2002, Wong et al. 1997
    # Km = 1.5 # For testing

    uptake_funs['EX_glc__D_e'] = lambda x: Vmax0 * x / (Km0 + x)

    # Glycerol:

    Vmax1 = 10*60 / 1000 * 56106 # umol/(min.mg[E]) --> mmol/(h.mmol[E])
                                # Voegele et al. 1992
                                # Pettigrew et al. 1987
    # Vmax1 = 10
    Km1 = 0.010 # mM, Voegele et al. 1992

    uptake_funs['EX_glyc_e'] = lambda x: Vmax1 * x / (Km1 + x)

    return uptake_funs

def prepare_model(in_model, S0_glyc, S0_glc):
    """
    Minimize glucose enzymes on glycerol uptake
    :param in_model:
    :return:
    """

    in_model.reactions.EX_glc__D_e.lower_bound = -1 * S0_glc
    in_model.reactions.EX_glyc_e.lower_bound = -1 * S0_glyc
    sol_glyc = in_model.optimize()

    # in_model.growth_reaction.lower_bound = sol_glyc.f * 0.9
    #
    # glycolysis_rxns = [x for x in in_model.reactions
    #                    if x.subsystem == 'Glycolysis/Gluconeogenesis' and
    #                    isinstance(x, EnzymaticReaction)]
    # glycolysis_enz = [item for rxn in glycolysis_rxns for item in rxn.enzymes]
    #
    # obj = sum([-1*x.concentration for x in glycolysis_enz])
    #
    # in_model.objective = obj
    #
    # sol_min = in_model.optimize()
    #
    # in_model.objective = in_model.growth_reaction
    # in_model.growth_reaction.lower_bound = 0
    #
    # return sol_min
    return sol_glyc


if __name__ == '__main__':

    timestep = 0.1

    S0_gly = 0 #mmol/L
    S1_gly = 0 #mmol/L

    S0_glc = 1e0 #mmol/L
    S1_glc = 10 #mmol/L
    # X0 = 0.1 #mmol[Cell]/L
    X0 = 0.05 #g[Cell]/L

    glc_step = lambda t,S,S0=S0_glc,S1=S1_glc: \
        S + S1 if abs(t - 1) <= timestep and S <= S0 else max(S, 0)
        # S1 if t > 1  else S0
    gly_step = lambda t,S,S0=S0_gly,S1=S1_gly: \
        S + S1 if abs(t - 1) <= timestep and S <= S0 else max(S, 0)

    glyc_glc_switch = {
        'EX_glc__D_e': glc_step,
        'EX_glyc_e'  : gly_step,
    }

    times = [x*timestep for x in range(1,10)]


    print("Test step fun:")
    print(times)
    for rxn, fun in glyc_glc_switch.items():
        values = [S0_glc]
        _ = [values.append(fun(t,values[-1])) for t in times]
        print(rxn, values)

    standard_solver_config(model)
    min_glycolysis_enz = prepare_model(model, S0_glyc=S0_gly, S0_glc=S0_glc)

    time_data = run_dynamic_etfl(model,
                                 timestep=timestep,
                                 tfinal=5,
                                 uptake_fun=get_uptake_funs(),
                                 medium_fun = glyc_glc_switch,
                                 S0={'EX_glyc_e':S0_gly, 'EX_glc__D_e':S0_glc},
                                 X0=X0,
                                 inplace=True,
                                 initial_solution = min_glycolysis_enz,
                                 chebyshev_include = [ForwardCatalyticConstraint,
                                                      BackwardCatalyticConstraint,
                                                      ExpressionCoupling]
                                 )

    plot_dynamics(model, time_data)
    time_data.to_csv('outputs/detfl_cheby.csv')
