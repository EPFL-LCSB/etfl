# -*- coding: utf-8 -*-

import pandas as pd
import bokeh.plotting as bp
from  bokeh.plotting import figure, show, output_file, curdoc

from collections import OrderedDict

from bokeh.palettes import Category10
from bokeh.io import export_svgs
from bokeh.layouts import row

from os.path import join as pjoin

bp.curdoc().clear()

data_dir = '../organism_data/info_ecoli'

from etfl.io.json import load_json_model
# scaling_model = load_json_model('models/RelaxedModel iJO1366_T1E1N1_431_enz_128_bins__20180926_124941.json')
scaling_model = load_json_model('models/iJO1366_EFL_431_enz_128_bins__20190121_080047.json')

def make_polygon(x, ymin, ymax):
    lx = x.tolist()
    lymin = ymin.tolist()
    lymax = ymax.tolist()

    return lx + lx[::-1], lymin + lymax[::-1]

def make_bars(x, ymin, ymax, width = None):
    data = pd.DataFrame([x,ymin,ymax]).T
    data.columns = ['x','bottom','top']
    data.dropna(how = 'any', inplace=True)

    if width is None:
        width = abs(x.diff()).min()/2

    top    = data['top']
    bottom = data['bottom']
    left   = data['x']-width/2
    right  = data['x']+width/2

    return [top, bottom, left, right]


def plot_growth_vs_uptake(fig, df, legend_text, color ='forestgreen'):
    fig.circle(df['uptake'], df['mu'], color=color, legend=legend_text)
    fig.line  (df['uptake'], df['mu'], color=color, legend=legend_text)

    # plot_var_vs_uptake(fig, df, 'mu', 1, legend_text, color)

def plot_var_vs_uptake(fig, df, var, sigma, legend_text, color, alpha = 0.5):

    # xs, ys = make_polygon(df['uptake'], df[var+'_lb'], df[var+'_ub'])
    # fig.patch(xs, ys, color=color, alpha=alpha, legend = legend_text)

    top, bottom, left, right = make_bars(df['uptake'], df[var+'_lb']*sigma, df[var+'_ub']*sigma, width = 0.5)
    fig.quad(top=top,
             bottom=bottom,
             left=left,
             right=right,
             color=color,
             alpha=alpha,
             legend=legend_text,
             line_width=2,
             line_alpha=2)


def plot_prop_vs_uptake(fig, df, prop, legend_text, color, plot_lines = True):
    fig.circle(df['uptake'], df[prop], color=color, legend=legend_text)
    if plot_lines:
        fig.line  (df['uptake'], df[prop], color=color, legend=legend_text)

def plot_prop_vs_growth(fig, df, prop, legend_text, color, plot_lines = True):
    fig.circle(df['mu'], df[prop], color=color, legend=legend_text)
    if plot_lines:
        fig.line  (df['mu'], df[prop], color=color, legend=legend_text)

# Neidhardt data
neidhardt_data = pd.read_excel(pjoin(data_dir,'neidhardt_tab2.xlsx'),
                               skiprows=range(0,6),
                               skipfooter=22)
mu_cols = ['mu=0.6','mu=1.0','mu=1.5','mu=2.0','mu=2.5']
neidhardt_data.columns = ['parameter','symbol','units'] + mu_cols +\
                          ['observed_parameters','footnote']
neidhardt_data.set_index('symbol', inplace=True)

Pc = neidhardt_data.loc['Pc (μg)'][mu_cols] # μg/10^9 cells
Rc = neidhardt_data.loc['Rc (μg)'][mu_cols] # μg/10^9 cells
Mc = neidhardt_data.loc['Mc (μg)'][mu_cols] # μg dry weight/10^9 cells

neidhardt_prel = (Pc/Mc).astype(float)
neidhardt_rrel = (Rc/Mc).astype(float)
neidhardt_mu = pd.Series(Pc.index.str.replace('mu=','')).astype(float)


def wrapper_plot_growth(data, color_offset=0):
    global mu_figure, e, model_name, df
    mu_figure = figure()
    for e, (model_name, df) in enumerate(data.items()):
        plot_growth_vs_uptake(mu_figure,
                              df,
                              model_name,
                              cmap[e+color_offset])
    mu_figure.legend.location = 'bottom_right'
    mu_figure.title.text = 'Growth rate vs glucose uptake'
    mu_figure.xaxis.axis_label = 'glucose uptake [mmol/(gDw.h)]'
    mu_figure.yaxis.axis_label = 'Growth rate [1/h]'
    mu_figure.output_backend = 'svg'

    return mu_figure

def wrapper_plot_var(data, var, color_offset=0):
    if var in ['EZ_rib', 'EZ_rnap']:
        fig = figure(y_axis_type='log')
    else:
        fig = figure()

    for e, (model_name, df) in enumerate(data.items()):

        if var.startswith('EZ_'):
            # it's an enzyme
            sigma = scaling_model.enzymes.get_by_id(var[3:]).scaling_factor
        elif var.startswith('MR_'):
            # it's an mRNA
            sigma = scaling_model.mrnas.get_by_id(var[3:]).scaling_factor
        else:
            sigma = 1

        plot_var_vs_uptake(fig,
                           df,
                           var,
                           sigma,
                           model_name,
                           cmap[e+color_offset],
                           alpha=0.8
                           )

    if var in ['EZ_dummy_enzyme']:
        fig.legend.location = 'top_right'
    elif var in ['EZ_rnap', 'EZ_rib']:
        fig.legend.location = 'bottom_right'
    else:
        fig.legend.location = 'top_left'

    fig.title.text = var
    fig.xaxis.axis_label = 'glucose uptake [mmol/(gDw.h)]'
    fig.yaxis.axis_label = var + ' concentration [μmol/gDw]'

    fig.output_backend = 'svg'

    return fig

def wrapper_plot_ratio(data, ratio, color_offset = 0):
    fig = figure()
    for e, (model_name, df) in enumerate(data.items()):
        plot_prop_vs_uptake(fig, df, ratio, model_name, cmap[e+color_offset])

    fig.title.text = ratio + ' vs glucose uptake'
    fig.xaxis.axis_label = 'glucose uptake [mmol/(gDw.h)]'
    fig.yaxis.axis_label = ratio + ' [g/gDw]'

    if 'mrna' in ratio:
        fig.legend.location = 'bottom_right'
    fig.output_backend = 'svg'

    return fig

def wrapper_plot_ratio_vs_growth(data,ratio,color_offset = 0):
    fig = figure()
    for e, (model_name, df) in enumerate(data.items()):
        plot_prop_vs_growth(fig, df, ratio, model_name, cmap[e], plot_lines=False)

    fig.title.text = ratio + ' vs growth rate'
    fig.xaxis.axis_label = 'Growth rate [1/h]'
    fig.yaxis.axis_label = ratio + ' [g/gDw]'

    if 'mrna' in ratio:
        fig.legend.location = 'bottom_right'
        fig.square(neidhardt_mu, neidhardt_rrel,
                   line_color=cmap[3],
                   fill_alpha=0,
                   size=10,
                   line_width=2,
                   legend='Neidhardt et al. measurements')
    elif 'prot' in ratio:
        fig.square(neidhardt_mu, neidhardt_prel,
                   line_color=cmap[3],
                   fill_alpha=0,
                   size=10,
                   line_width=2,
                   legend='Neidhardt et al. measurements')
    fig.output_backend = 'svg'

    return fig

if __name__ == '__main__':

    cmap = Category10[10]

    model_data1 = OrderedDict({
        'EFL': pd.read_csv('outputs/benchmark_EFL_25.csv'),
        'ETFL': pd.read_csv('outputs/benchmark_ETFL_25.csv'),
        'vEFL': pd.read_csv('outputs/benchmark_vEFL_25.csv'),
        'vETFL': pd.read_csv('outputs/benchmark_vETFL_25.csv'),
    })
    model_data2 = OrderedDict({
        'vETFL': model_data1['vETFL'],
        'vETFL65': pd.read_csv('outputs/benchmark_vETFL65_25.csv'),
        'vETFL, inferred enzymes': pd.read_csv('outputs/benchmark_vETFL_infer_25.csv'),
        'vETFL65, inferred enzymes': pd.read_csv('outputs/benchmark_vETFL65_infer_25.csv'),
    })
    model_data_old = OrderedDict({
        # 'T0E1N0': pd.read_csv('outputs/benchmark_T0E1N0.csv'),
        # 'T0E1N1': pd.read_csv('outputs/benchmark_T0E1N1.csv'),
        # 'T1E1N0': pd.read_csv('outputs/benchmark_T1E1N0.csv'),
        # 'T1E1N1': pd.read_csv('outputs/benchmark_T1E1N1.csv'),
        # 'T0E1N1_inferred_mean': pd.read_csv('outputs/benchmark_T0E1N1_inferred_mean.csv'),
        # 'T0E1N1_inferred_median': pd.read_csv('outputs/benchmark_T0E1N1_inferred_median.csv'),
    })


    OFFSET = len(model_data1)-1

    # Plot growth rate vs glucose uptake


    mu_figure1 = wrapper_plot_growth(data = model_data1)
    mu_figure2 = wrapper_plot_growth(data = model_data2,
                                     color_offset=OFFSET)
    growth_fig = row([mu_figure1,mu_figure2])
    filename  = 'plots/benchmark_growth'
    output_file(filename + '.html')
    show(growth_fig)
    export_svgs(growth_fig, filename=filename + '.svg')

    curdoc().clear()

    # Plot variability
    variables_to_plot = ['EZ_rib',
                         'EZ_rnap',
                         # 'EZ_dummy_enzyme',
                         # 'MR_dummy_gene',
                         ]

    for var in variables_to_plot:
        var_fig1 = wrapper_plot_var(model_data1,var)
        var_fig2 = wrapper_plot_var(model_data2,var,color_offset=OFFSET)
        fig = row([var_fig1,var_fig2])

        filename = 'plots/benchmark_{}'.format(var)
        output_file(filename + '.html')
        show(fig)
        export_svgs(fig, filename=filename + '.svg')
        curdoc().clear()

    # mrna / ribosome content

    for ratio in ['mrna_ratio','prot_ratio']:
        ratio_fig1 = wrapper_plot_ratio(model_data1, ratio)
        ratio_fig2 = wrapper_plot_ratio(model_data2, ratio, color_offset=OFFSET)
        fig = row([ratio_fig1, ratio_fig2])

        filename = 'plots/benchmark_{}'.format(ratio)
        output_file(filename + '.html')
        show(fig)
        export_svgs(fig, filename=filename + '.svg')
        curdoc().clear()


    # mrna / ribosome content

    for ratio in ['mrna_ratio','prot_ratio']:
        ratio_fig3 = wrapper_plot_ratio_vs_growth(model_data1, ratio)
        ratio_fig4 = wrapper_plot_ratio_vs_growth(model_data1, ratio, color_offset=OFFSET)

        fig = row([ratio_fig3,ratio_fig4])

        filename = 'plots/benchmark_{}_vs_mu'.format(ratio)
        output_file(filename + '.html'.format(ratio))

        show(fig)
        export_svgs(fig, filename=filename + '.svg')
        curdoc().clear()
