from etfl.tests.small_model import create_etfl_model
from etfl.analysis.dynamic import run_dynamic_etfl
from etfl.optim.config import standard_solver_config
from etfl.core.reactions import EnzymaticReaction
from etfl.optim.constraints import ForwardCatalyticConstraint, \
    BackwardCatalyticConstraint, ExpressionCoupling

import bokeh.plotting as bp
from bokeh.palettes import Category10
from bokeh.layouts import row, gridplot
from bokeh.models import LinearAxis, Range1d

from bokeh.io import export_svgs

from math import exp
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

model = load_json_model('models/SlackModel iJO1366_vETFL_440_enz_128_bins__20190227_141232.json')
# model = load_json_model('models/SlackModel iJO1366_vETFL_2084_enz_128_bins__20190122_170118.json')
# model = load_json_model('models/iJO1366_vEFL_431_enz_128_bins__20190121_090316.json')
# model = load_json_model('models/SlackModel iJO1366_vETFL_431_enz_128_bins__20190122_145700.json')
# model = load_json_model('models/iJO1366_T0E1N1_431_enz_128_bins__20180926_135704.json')
# model = load_json_model('models/RelaxedModel iJO1366_T1E1N1_431_enz_128_bins__20180926_124941.json')

model_tag = 'dvETFL_deg_glc_lcts_alldelta_noexpc'

# for rxn in model.translation_reactions:
#     rxn.subsystem = 'Translation'
# for rxn in model.transcription_reactions:
#     rxn.subsystem = 'Transcription'


exchanges = [
    'EX_glc__D_e',
    'EX_ac_e',
    'EX_o2_e',
    'EX_lcts_e']

concentrations = ['X'] + ['S_'+x for x in exchanges]

flux_ubs = [('act_'+x,'ub_'+x) for x in exchanges]

def get_enzymes_of_subsystem(model, subsystem):
    reactions = [x for x in model.reactions if subsystem.lower() in x.subsystem.lower()]

    enzymes = [item for rxn in reactions if hasattr(rxn,'enzymes') and rxn.enzymes is not None
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

def plot_subsystems(p,t,model,time_data):

    subsystems = list(set((x.subsystem for x in model.reactions
                           if x.subsystem is not None and x.subsystem)))

    all_enz = dict()
    for sub in subsystems:
        these_enzymes = get_enzymes_of_subsystem(model, sub)
        all_enz[sub] = get_prot_total(time_data, these_enzymes)

    data = pd.DataFrame.from_dict(all_enz, orient = 'columns')
    data = data.reindex(data.max().sort_values(ascending=False).index, axis=1)
    chosen_subs = data.columns[:10]
    data = data[chosen_subs]

    data['t'] = t
    # data.index = t
    # data = data.reset_index()
    from bokeh.palettes import Category10
    from bokeh.core.properties import value
    colors = Category10[len(data.columns)-1]

    last_y = 0*t
    times = list(t) + list(t[::-1])

    for e,sub in enumerate(chosen_subs):
        this_y = data[sub] + last_y
        ys = list(this_y) + list(last_y[::-1])

        p.patch(x = times, y = ys, color = colors[e], legend = sub)

        last_y = this_y

    # p.vbar_stack(chosen_subs,
    #              x='t',
    #              width=t.diff().min(),
    #              color=colors,
    #              source=data,
    #              legend=[value(x) for x in chosen_subs],
    #              )
    p.legend.location = 'top_left'

    return data

def plot_growth(p,t,time_data):

    plot_var(p,t,time_data, 'X', color='black')
    p.yaxis.axis_label = "Cell concentration $[g.L^{-1}]$"
    p.y_range.start = 0.9 * time_data.loc['X'].min()
    p.y_range.end =  1.09 * time_data.loc['X'].max()

    # Setting the second y axis range name and range
    p.extra_y_ranges = {"mu": Range1d()}

    # Adding the second axis to the plot.
    p.add_layout(LinearAxis(y_range_name="mu",
                            axis_label='growth rate $[h^{-1}]$'), 'right')

    plot_var(p,t,time_data, 'mu', color='gray', y_range_name='mu')

    return p


def plot_dynamics(model, time_data):

    bp.curdoc().clear()

    # time_data += model.solver.configuration.tolerances.feasibility

    p1  = bp.figure()#y_axis_type = 'log')
    p2  = bp.figure()#y_axis_type = 'log')
    p3  = bp.figure(y_axis_type = 'log')
    p4  = bp.figure(y_axis_type = 'log')
    ps  = bp.figure()
    px  = bp.figure()

    t = time_data.loc['t']

    vars_to_obs = ['MR_b4025', 'MR_b1779', 'EZ_rib', 'EZ_rnap']

    cmap = Category10[10]

    for e, the_var in enumerate(concentrations):
        plot_var(p1, t, time_data, the_var, cmap[e])

    for e, (the_flux, the_ub) in enumerate(flux_ubs):
        plot_ub (p2, t, time_data, the_ub  , color=cmap[e])
        plot_var(p2, t, time_data, the_flux, color=cmap[e])

    p2.legend.location = 'bottom_right'
    plot_var(p2,t,time_data, 'mu', color=cmap[e+1])

    plot_growth(px,t,time_data)

    glycolysis_enzymes = get_enzymes_of_subsystem(model, 'Glycolysis')
    ppp_enzymes = get_enzymes_of_subsystem(model, 'pentose phosphate')
    tca_enzymes = get_enzymes_of_subsystem(model, 'Citric Acid Cycle')


    for e, the_var in enumerate(vars_to_obs):
        plot_var(p3, t, time_data, the_var, cmap[e])

    time_data.loc['glycolysis',:] = get_prot_total(time_data, glycolysis_enzymes)
    time_data.loc['ppp',:] = get_prot_total(time_data, ppp_enzymes)
    time_data.loc['tca',:] = get_prot_total(time_data, tca_enzymes)
    plot_var(p4, t, time_data, 'glycolysis', cmap[e + 1])
    plot_var(p4, t, time_data, 'ppp', cmap[e + 2])
    plot_var(p4, t, time_data, 'tca', cmap[e + 3])

    effluxes = time_data.loc[ time_data.index.str.startswith('EX_')
                              & time_data.index.str.endswith('_e') ]

    effluxes['score'] = effluxes.std(axis=1) / effluxes.mean(axis=1)

    effluxes = effluxes.reindex(
        effluxes['score'].sort_values(ascending = False).index,
        axis=0)

    # for e,(the_var, values) in enumerate(effluxes[:5].T.iteritems()):
    #
    #     pex.line(t, values,
    #              color = cmap[e],
    #              line_width = 2,
    #              legend=the_var)

    plot_subsystems(ps,t,model,time_data)

    for p in [p1,p2,p3,p4,ps,px]:
        #     p.legend.location = 'top_left'
        p.output_backend = 'svg'

    gp = gridplot([[p1,p3],[p2,p4],[ps,px]])
    bp.output_file('plots/{}.html'.format(model_tag))
    bp.show(gp)
    export_svgs(gp, filename='plots/{}.svg'.format(model_tag))

def plot_var(p, t, time_data, the_var, color, **kwargs):

    p.circle(t, time_data.loc[the_var],
             line_color=color,
             fill_alpha=0,
             legend=the_var,
             size=10,
             line_width=2,
             **kwargs
             )

    p.line(t, time_data.loc[the_var],
           line_color=color,
           line_width=2,
           # legend=the_var,
           **kwargs
           )

def plot_ub(p, t, time_data, the_var, color):

    p.square(t, time_data.loc[the_var],
              # marker='x',
              line_color=color,
              fill_alpha=0,
              # legend=the_var,
              size=10,
              line_width=2
              )

def get_uptake_funs():

    uptake_funs = dict()

    # Glucose:

    # Vmax0 = 10 # mmol/(h.mmol[E]) Mahadevan et al. 2002
    Vmax0 = 1
    Km0 = 0.015 # mM, Mahadevan et al. 2002, Wong et al. 1997

    uptake_funs['EX_glc__D_e'] = lambda x: Vmax0 * x / (Km0 + x)

    # Lactose:
    # Olsen, S. G., and R. J. Brooker.
    # "Analysis of the structural specificity of the lactose permease toward sugars."
    # Journal of Biological Chemistry 264.27 (1989): 15982-15987.
    # http://www.jbc.org/content/264/27/15982.full.pdf+html
    # TODO: Stop assuming 1mgProt/gDW
    # Vmax_lac = 210 /1000 * 60 *1  # nmol/min/mgProt * mmol/nmol * min/h * mgProt/gDW
    Vmax_lac = 1
    Km_lac = 1.3 # mM

    uptake_funs['EX_lcts_e'] = lambda x: Vmax_lac * x / (Km_lac + x)

    # O2
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC134846/
    # Alexeeva S, Hellingwerf KJ, Teixeira de Mattos MJ.
    # Quantitative assessment of oxygen availability: perceived aerobiosis and
    # its effect on flux distribution in the respiratory chain of
    # Escherichia coli.
    # J Bacteriol. 2002;184(5):1402-6.
    # The cytochrome bo oxidase has a Km for O2 of 2 × 10−4 mM and a Vmax of
    # 6.6 × 10−2 mmol of O2·nmol of cytochrome o−1·h−1 (18). This corresponds,
    # at the measured rDOT value of 1.6 × 10−2 mM, to a cytochrome bo oxidase
    # content of 73 nmol of protein·g (dry weight)−1
    # VmaxO2 = 6.6*1e-2 *73 # mmol/h
    VmaxO2 = 15 # mmol/h
    KmO2 = 2 * 1e-4 #mM
    uptake_funs['EX_o2_e'] = lambda x: VmaxO2 * x /(KmO2 + x)

    # Acetate
    uptake_funs['EX_ac_e'] = lambda x: 15 if x>0 else 0

    return uptake_funs

def get_uptake_enzymes(model):
    pass


def prepare_model(in_model, S0, uptake_fun):
    """
    Minimize glucose enzymes on glycerol uptake
    :param in_model:
    :return:
    """

    # for uptake_flux, kinfun in uptake_fun.items():
    #     in_model.reactions.get_by_id(uptake_flux).lower_bound = \
            # -1 * kinfun(S0[uptake_flux])
    model.reactions.EX_glc__D_e.lower_bound = -8
    model.reactions.EX_lcts_e.lower_bound = -0#8
    model.reactions.EX_o2_e.lower_bound = -5
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
    # timestep = 0.05
    # timestep = 0.02
    epsilon = timestep/100

    S0_o2 = 0  # mmol/L
    S1_o2 = 0.21  # mmol/L
    kla_o2 = 7.5  # h^-1

    S0_glc = 0.5 #mmol/L
    S1_glc = 10 #mmol/L

    S0_lac = 1 #mmol/L
    S1_lac = 10 #mmol/L

    S0_ac = 0


    # X0 = 0.1 #mmol[Cell]/L
    X0 = 0.05 #g[Cell]/L

    glc_fun = lambda t, S, S0=S0_glc, S1=S1_glc: \
        max(S, 0)
        # S + S1 if abs(t - 1) <= timestep+epsilon and S <= S0 else max(S, 0)
    # S1 if t > 1  else S0


    lac_fun = lambda t, S, S0=S0_lac, S1=S1_lac: \
        max(S, 0)

    o2_fun = lambda t, S, S0=S0_o2, S1=S1_o2: \
        S1 # if t > 1  else S0
    # S + S1 if abs(t - 1) <= timestep+epsilon and S <= S0 else max(S, 0)

    glc_free = lambda t,S : max(S0_glc,0)
    # o2_diff = lambda t,S,S0=S0_o2,kla=kla_o2 : max(S + kla*(S0-S),0)
    o2_diff = lambda t,S,S0=S1_o2,kla=kla_o2 : max(S0 - (S0-S)*exp(-kla*timestep),0)

    ac_fun = lambda t,S : max(S,0)

    conc_funs = {
        # 'EX_glc__D_e': glc_step,
        # 'EX_o2_e'    : o2_step,
        'EX_glc__D_e': glc_fun,
        'EX_lcts_e': lac_fun,
        'EX_o2_e'    : o2_diff,
        'EX_ac_e'    : ac_fun,
    }

    S0 = {'EX_o2_e': S1_o2,
          'EX_glc__D_e': S0_glc,
          'EX_lcts_e': S0_lac,
          'EX_ac_e': S0_ac}

    uptake_enz = {
        'EX_glc__D_e': ['GLCDpp', 'GLCabcpp','GLCptspp', 'GLCt2pp'],
        'EX_lcts_e'  : ['LACZpp', 'LCTSt3ipp', 'LCTStpp'],
        'EX_ac_e'    : ['ACt4pp']
        }

    times = [x*timestep for x in range(1,10)]


    print("Test step fun:")
    print(times)
    for rxn, fun in conc_funs.items():
        values = [S0[rxn]]
        _ = [values.append(fun(t,values[-1])) for t in times]
        print(rxn, values)

    standard_solver_config(model)
    model.solver.configuration.verbosity = 0

    ini_sol = prepare_model(model, S0=S0, uptake_fun = get_uptake_funs())

    time_data = run_dynamic_etfl(model,
                                 timestep=timestep,
                                 tfinal=3,
                                 uptake_fun=get_uptake_funs(),
                                 uptake_enz=uptake_enz,
                                 medium_fun =conc_funs,
                                 S0=S0,
                                 X0=X0,
                                 inplace=True,
                                 initial_solution = ini_sol,
                                 chebyshev_include = [ForwardCatalyticConstraint,
                                                      BackwardCatalyticConstraint,
                                                      ExpressionCoupling,
                                                      ]
                                 )

    plot_dynamics(model, time_data)
    time_data.to_csv('outputs/{}.csv'.format(model_tag))
