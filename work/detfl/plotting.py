import pandas as pd

import bokeh.plotting as bp
from bokeh.layouts import column, gridplot
from bokeh.palettes import Category10, magma, plasma
from bokeh.models import Legend, Range1d, LinearAxis

from os import makedirs
from os.path import join, exists

from pytfa.utils.logger import get_timestr

AXIS_FONT_SIZE = "25pt"
LEGEND_FONT_SIZE = "25pt"
LINE_WIDTH = 6

colorblind4 = ['#D81B60',
               '#1E88E5',
               '#FFC107',
               '#004D40',
               ]

def get_mrna_total(time_data, mrnas):
    return time_data.loc[[x.variable.name * x.scaling_factor
                          for x in mrnas]].sum()

def get_prot_total(time_data, enzymes, mass_mode = True):
    enz = time_data.loc[[x.variable.name for x in enzymes]]
    if mass_mode:
        sigma=1
    else:
        sigma =  pd.Series([x.scaling_factor for x in enzymes], index = enz.index)

    return (sigma*enz.T).T.sum()


def summarize_model(model,time_data,groups,
                    output_path='.', model_tag='', backend ='png'):

    summary_plots = dict()
    detailed_plots = dict()
    cleaned_groups = dict()

    # Cleanup missing vars
    for key, varnames in groups.items():
        cleaned_groups[key] = [x for x in varnames if x in time_data.index or x == 'total']

    groups = cleaned_groups
    total_groups = {k:v for k,v in groups.items() if 'total' in v and 'pathway_enz' in k}

    for key,data_type in groups.items():
        summary_plots[key] = make_summary_plots(model,time_data,key,data_type)

        detailed_plots[key] = make_detailed_plots(model,time_data,data_type,backend)

    summary_plots['growth'] = make_growth_plot(model,time_data)

    summary_plots['mass_ratios'] = plot_mass(model,time_data)

    if model is not None:
        summary_plots['subsystems'] = plot_subsystems(model,time_data)

    summary_plots['total'] = make_total_plot(model, time_data, total_groups)

    output_folder = join(output_path)#,model_tag)#+'_'+get_timestr())

    if not exists(output_folder):
        makedirs(output_folder)

    bp.curdoc().clear()
    bp.output_file(join(output_folder,'summary.html'), title = model_tag)

    for p in summary_plots.values():
        p.output_backend = backend

    bp.show(column(list(summary_plots.values())))

    for key,this_dp in detailed_plots.items():
        bp.curdoc().clear()
        bp.output_file(join(output_folder,'{}.html'.format(key)))
        for p in this_dp.children:
            try:
                p.output_backend = backend
            except AttributeError:
                # Fails for gridplots
                pass
                # from bokeh.layouts import GridBox
                # if isinstance(p,GridBox):
                #     p.children[0][0].output_backend = backend
                # else: # Toolbox bar etc..
                #     pass
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

def make_total_plot(model, time_data, total_groups):

    t = time_data.loc['t']
    tot_y = pd.DataFrame(columns = t.index)

    for e,(groupname,vars) in enumerate(total_groups.items()):
        ys = pd.DataFrame.from_dict(
            {the_var:get_y(the_var, time_data) for the_var in vars},
            orient='index')
        tot_y.loc[groupname] = ys.sum()

    p = plot_lines(t, tot_y, title='')
    return p




def make_detailed_plots(model, time_data, data_type, backend):

    t = time_data.loc['t']

    p = []
    for the_var in data_type:
        #is there a reverse var ?
        y = get_y(the_var, time_data)

        this_p = plot_line(t,y,enhance_varnames(the_var))
        p.append(this_p)

        this_p.xaxis.axis_label = 'time [h]'
        if the_var.startswith('EZ_'):
            this_p.yaxis.axis_label = 'Enzyme mass [g/gDW]'

            # return column(p)
    for this_p in p:
        this_p.output_backend = backend

    return gridplot(p, ncols=3)

def make_growth_plot(model,time_data):
    t = time_data.loc['t']
    y1 = time_data.loc['X']
    y2 = time_data.loc['mu']

    p = bp.figure(width=1000)
    p.title.text = 'Growth, Cell concentration over time'
    p.line(t,y1, color='black', line_width=LINE_WIDTH)

    p.y_range.start = -0.05 * y1.min()
    p.y_range.end   =  1.05 * y1.max()

    # Setting the second y axis range name and range
    p.extra_y_ranges = {"mu": Range1d(start=-0.05*y2.min(), end = 1.05*y2.max())}

    # Adding the second axis to the plot.
    p.add_layout(LinearAxis(y_range_name="mu",
                            axis_label='growth rate $[h^{-1}]$'), 'right')

    p.line(t,y2, color='grey',  line_width=LINE_WIDTH,
           line_dash = 'dashed', y_range_name='mu')

    p.xaxis.major_label_text_font_size = AXIS_FONT_SIZE
    p.yaxis.major_label_text_font_size = AXIS_FONT_SIZE
    # p.legend.label_text_font_size = LEGEND_FONT_SIZE

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

def enhance_varnames(s):
    """ Removes pieces of varnames to make them more legible"""
    return s.replace('EZ_','').replace('_MONOMER',''.split('_mod_')[0])

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
        the_legend = enhance_varnames(labels[e])
        c = p.line(t,y, color=colors[e], line_width=LINE_WIDTH)#, legend=labels[e])
        legend_it.append((the_legend, [c]))

    if total:
        tot_y = ys.sum()
        c = p.line(t, tot_y, color='grey', line_width=LINE_WIDTH, line_dash='dashed')
        legend_it.append(('total', [c]))


    legend = Legend(items=legend_it, location=(0, 0))

    p.add_layout(legend, 'right')
    p.xaxis.major_label_text_font_size = AXIS_FONT_SIZE
    p.yaxis.major_label_text_font_size = AXIS_FONT_SIZE
    # p.legend.label_text_font_size = LEGEND_FONT_SIZE

    return p


def plot_subsystems(model,time_data,compact=True):
    t = time_data.loc['t']
    if compact:
        p = bp.figure(width=1000)
    else:
        p = bp.figure(width=1000, height = 800)

    subsystems = list(set((x.subsystem for x in model.reactions
                           if x.subsystem is not None and x.subsystem)))

    all_enz = dict()
    for sub in subsystems:
        these_enzymes = get_enzymes_of_subsystem(model, sub)
        if len(these_enzymes) ==0:
            continue
        all_enz[sub] = get_prot_total(time_data, these_enzymes)

    data = pd.DataFrame.from_dict(all_enz, orient = 'columns')
    data = data.reindex(data.max().sort_values(ascending=False).index, axis=1)
    data.to_csv('tmp.csv')

    if compact:
        chosen_subs = data.columns[:10]
        colors = Category10[len(chosen_subs)]
    else:
        chosen_subs = data.columns
        colors = plasma(len(chosen_subs))

    data = data[chosen_subs]

    data['t'] = t

    last_y = 0*t
    times = list(t) + list(t[::-1])
    legend_it = []
    for e,sub in enumerate(chosen_subs):
        this_y = data[sub] + last_y
        ys = list(this_y) + list(last_y[::-1])

        c = p.patch(x = times, y = ys, color = colors[e])#, legend = sub)
        legend_it.append((sub, [c]))


        last_y = this_y

    p.legend.location = 'top_left'
    legend = Legend(items=legend_it[::-1], location=(0, 0))

    p.add_layout(legend, 'right')

    return p

def plot_mass(model,time_data):

    # p = bp.figure(width=1000)

    t = time_data.loc['t']
    # prot_mass = time_data.loc['IV_prot_ggdw']
    # mrna_mass = time_data.loc['IV_mrna_ggdw']
    # dna_mass  = time_data.loc['IV_dna_ggdw']
    ys = time_data.loc[['IV_prot_ggdw','IV_mrna_ggdw','IV_dna_ggdw']]

    # colors = Category10[3]
    #
    # last_y = 0*t
    # times = list(t) + list(t[::-1])
    # legend_it = []

    p = plot_lines(t, ys, title = 'Mass ratios',
                   total=True)

    # for e,data in enumerate([prot_mass,mrna_mass,dna_mass]):
    #     this_y = data + last_y
    #     ys = list(this_y) + list(last_y[::-1])
    #     c = p.patch(x = times, y = ys, color = colors[e])#, legend = sub)
    #     legend_it.append((data.name, [c]))
    #     last_y = this_y
    #
    # p.legend.location = 'top_left'
    # legend = Legend(items=legend_it[::-1], location=(0, 0))
    #
    # p.add_layout(legend, 'right')

    return p


def get_enzymes_of_subsystem(model, subsystem):
    reactions = [x for x in model.reactions if subsystem.lower() in x.subsystem.lower()]

    enzymes = [item for rxn in reactions if hasattr(rxn,'enzymes') and rxn.enzymes is not None
               for item in rxn.enzymes
               ]

    return enzymes

if __name__ == '__main__':

    if not exists('plots'):
        makedirs('plots/tmp')

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
    summarize_model(model, time_data, groups, output_path='plots/tmp',
                    model_tag=model_tag)
