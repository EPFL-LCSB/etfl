import pandas as pd

import bokeh.plotting as bp
from bokeh.layouts import column
from bokeh.palettes import Category10, magma
from bokeh.models import Legend, Range1d, LinearAxis

from os import makedirs
from os.path import join, exists

from pytfa.utils.logger import get_timestr



def get_mrna_total(time_data, mrnas):
    return time_data.loc[[x.variable.name * x.scaling_factor
                          for x in mrnas]].sum()

def get_prot_total(time_data, enzymes):
    enz = time_data.loc[[x.variable.name for x in enzymes]]
    sigma =  pd.Series([x.scaling_factor for x in enzymes], index = enz.index)

    return (sigma*enz.T).T.sum()


def summarize_model(model,time_data,groups,
                    output_path='.', model_tag=''):

    summary_plots = dict()
    detailed_plots = dict()
    cleaned_groups = dict()

    # Cleanup missing vars
    for key, varnames in groups.items():
        cleaned_groups[key] = [x for x in varnames if x in time_data.index or x == 'total']

    groups = cleaned_groups

    for key,data_type in groups.items():
        summary_plots[key] = make_summary_plots(model,time_data,key,data_type)

        detailed_plots[key] = make_detailed_plots(model,time_data,data_type)

    summary_plots['growth'] = make_growth_plot(model,time_data)

    output_folder = join(output_path,model_tag)#+'_'+get_timestr())

    if not exists(output_folder):
        makedirs(output_folder)

    bp.curdoc().clear()
    bp.output_file(join(output_folder,'summary.html'))
    bp.show(column(list(summary_plots.values())))

    for key,this_dp in detailed_plots.items():
        bp.curdoc().clear()
        bp.output_file(join(output_folder,'{}.html'.format(key)))
        bp.show(this_dp)

def make_summary_plots(model, time_data, key, data_type, total=False):

    t = time_data.loc['t']

    try:
        total = data_type.pop(data_type.index('total'))
    except ValueError:
        total = False

    ys = pd.DataFrame.from_dict(
        {the_var:get_y(the_var, time_data) for the_var in data_type},
        orient='index')
    p = plot_lines(t,ys,key, total)

    return p

def make_detailed_plots(model, time_data, data_type):

    t = time_data.loc['t']

    p = []
    for the_var in data_type:
        #is there a reverse var ?
        y = get_y(the_var, time_data)

        p.append(plot_line(t,y,the_var))

    return column(p)

def make_growth_plot(model,time_data):
    t = time_data.loc['t']
    y1 = time_data.loc['X']
    y2 = time_data.loc['mu']

    p = bp.figure(width=1000)
    # p.title.text = 'Growth, Cell concentration over time'
    p.line(t,y1, color='black', line_width=2)

    p.y_range.start = -0.05 * y1.min()
    p.y_range.end   =  1.05 * y1.max()

    # Setting the second y axis range name and range
    p.extra_y_ranges = {"mu": Range1d(start=-0.05*y2.min(), end = 1.05*y2.max())}

    # Adding the second axis to the plot.
    p.add_layout(LinearAxis(y_range_name="mu",
                            axis_label='growth rate $[h^{-1}]$'), 'right')

    p.line(t,y2, color='grey',  line_width=2, line_dash = 'dashed', y_range_name='mu')

    return p

def get_y(the_var, time_data):
    reverse_var = time_data.index.str.contains(the_var + '_reverse')
    if sum(reverse_var) > 0:
        y = time_data.loc[the_var] - time_data[reverse_var].iloc[0]
    else:
        y = time_data.loc[the_var]
    return y


def plot_line(t,y,label, color = 'black'):

    low  = min(-0.05*max(y), 1.05*min(y))
    high = max(-0.05*min(y), 1.05*max(y))

    if low == high:
        low  = -1
        high =  1

    p = bp.figure(height=300, width=300, y_range=(low,high))
    p.title.text = label
    p.line(t,y, color=color, line_width=2)

    return p


def plot_lines(t,ys,title, total = False):

    p = bp.figure(width = 1000)#height=300, width=300)
    p.title.text = title
    labels = list(ys.index)
    if len(ys) <= 10:
        colors = Category10[10]
    else:
        colors = magma[len(ys)]

    legend_it = []

    for e,(row,y) in enumerate(ys.iterrows()):
        c = p.line(t,y, color=colors[e], line_width=2)#, legend=labels[e])
        legend_it.append((labels[e], [c]))

    if total:
        tot_y = ys.sum()
        c = p.line(t, tot_y, color='grey', line_width=2, line_dash='dashed')
        legend_it.append(('total', [c]))


    legend = Legend(items=legend_it, location=(0, 0))

    p.add_layout(legend, 'right')

    return p


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
    colors = Category10[len(data.columns)-1]

    last_y = 0*t
    times = list(t) + list(t[::-1])

    for e,sub in enumerate(chosen_subs):
        this_y = data[sub] + last_y
        ys = list(this_y) + list(last_y[::-1])

        p.patch(x = times, y = ys, color = colors[e], legend = sub)

        last_y = this_y

    p.legend.location = 'top_left'

    return data

if __name__ == '__main__':

    if not exists('plots'):
        makedirs('plots')

    # time_data_path = 'data/detfl_cheby_monod_vETFL_vmax_mixed_expc_lcts2_with_degradation.csv'
    time_data_path = 'data/detfl_cheby_monod_vETFL_vmax_mixed_ini_lcts2_with_degradation.csv'
    # time_data_path = 'outputs/dvETFL_deg_glc_lcts_alldelta_noexpc.csv'

    # model_tag = 'glc_lcts_no_expc_deg'
    model_tag = 'glc_lcts_mixed_no_expc_deg'

    time_data = pd.read_csv(time_data_path, header=0, index_col=0)

    model = None

    fluxes = ['EX_glc__D_e','EX_lcts_e','EX_o2_e','EX_ac_e']
    glc_fluxes = ['G1PPpp', 'GLCDpp', 'GLCabcpp', 'GLCptspp', 'GLCt2pp', 'TREHpp']
    lcts_fluxes = ['LACZpp', 'LCTSt3ipp', 'LCTStpp']

    lcts_enzymes = ['EZ_'+ x for x in ['LACZpp_EG12013_MONOMER',
                                       'LCTSt3ipp_B2170_MONOMER',
                                       'LCTSt3ipp_B0070_MONOMER',
                                       'LCTSt3ipp_YDEA_MONOMER',
                                       'LCTStpp_LACY_MONOMER',
                                      ]]
    glc_enzymes = ['EZ_'+ x for x in ['G1PPpp_GLUCOSE_1_PHOSPHAT_CPLX',
                                      'GLCDpp_GLUCDEHYDROG_MONOMER_mod_pqq',
                                      'GLCDpp_G6437_MONOMER_mod_ca2_mod_pqq',
                                      'GLCabcpp_ABC_18_CPLX',
                                      'GLCptspp_CPLX_164',
                                      'GLCptspp_CPLX_157',
                                      'GLCt2pp_GALP_MONOMER',
                                      'TREHpp_TREHALAPERI_MONOMER']]

    species = ['S_EX_' + x for x in ['glc__D_e','lcts_e','o2_e','ac_e']]
    groups = {'fluxes': fluxes,
              'glc_enzymes': glc_enzymes,
              'lcts_enzymes': lcts_enzymes,
              'glc_fluxes': glc_fluxes,
              'lcts_fluxes': lcts_fluxes,
              'species': species}
    summarize_model(model, time_data, groups, output_path='plots',
                    model_tag=model_tag)
