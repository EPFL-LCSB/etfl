# -*- coding: utf-8 -*-

import pandas as pd
import bokeh.plotting as bp
from  bokeh.plotting import figure, show, output_file, curdoc

from collections import OrderedDict

from bokeh.palettes import Category10

from os.path import join as pjoin

bp.curdoc().clear()

data_dir = '../organism_data/info_ecoli'

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

    plot_var_vs_uptake(fig, df, 'mu', legend_text, color)

def plot_var_vs_uptake(fig, df, var, legend_text, color, alpha = 0.5):

    # xs, ys = make_polygon(df['uptake'], df[var+'_lb'], df[var+'_ub'])
    # fig.patch(xs, ys, color=color, alpha=alpha, legend = legend_text)

    top, bottom, left, right = make_bars(df['uptake'], df[var+'_lb'], df[var+'_ub'], width = 0.5)
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

if __name__ == '__main__':

    cmap = Category10[10]

    model_data = OrderedDict({
        'T0E1N0': pd.read_csv('outputs/benchmark_T0E1N0.csv'),
        'T0E1N1': pd.read_csv('outputs/benchmark_T0E1N1.csv'),
        'T1E1N0': pd.read_csv('outputs/benchmark_T1E1N0.csv'),
        'T1E1N1': pd.read_csv('outputs/benchmark_T1E1N1.csv'),
    })


    # Plot growth rate vs glucose uptake

    mu_figure = figure()

    for e, (model_name, df) in enumerate(model_data.items()):
        plot_growth_vs_uptake(mu_figure,
                              df,
                              model_name,
                              cmap[e])

    mu_figure.legend.location = 'bottom_right'
    mu_figure.title.text = 'Growth rate vs glucose uptake'
    mu_figure.xaxis.axis_label = 'glucose uptake [mmol/(gDw.h)]'
    mu_figure.yaxis.axis_label = 'Growth rate [1/h]'

    output_file('plots/benchmark_growth.html')
    show(mu_figure)
    curdoc().clear()

    # Plot variability
    variables_to_plot = ['EZ_rib',
                         'EZ_rnap',
                         # 'EZ_dummy_enzyme',
                         # 'MR_dummy_gene',
                         ]
    figures = {}

    for var in variables_to_plot:
        if var in ['EZ_rib', 'EZ_rnap']:
            fig = figure(y_axis_type='log')
        else:
            fig = figure()

        for e, (model_name, df) in enumerate(model_data.items()):
            plot_var_vs_uptake( fig,
                                df,
                                var,
                                model_name,
                                cmap[e],
                                alpha=0.8
                                )

        output_file('plots/benchmark_{}.html'.format(var))

        if var in ['EZ_dummy_enzyme']:
            fig.legend.location = 'top_right'
        elif var in ['EZ_rnap', 'EZ_rib']:
            fig.legend.location = 'bottom_right'
        else:
            fig.legend.location = 'top_left'

        fig.title.text = var
        fig.xaxis.axis_label = 'glucose uptake [mmol/(gDw.h)]'
        fig.yaxis.axis_label = var + ' concentration [μmol/gDw]'



        show(fig)
        curdoc().clear()
        figures[model_name] = fig

    # mrna / ribosome content

    for ratio in ['mrna_ratio','prot_ratio']:
        fig = figure()
        for e, (model_name, df) in enumerate(model_data.items()):
            plot_prop_vs_uptake(fig,df,ratio,model_name,cmap[e])

        fig.title.text = ratio + ' vs glucose uptake'
        fig.xaxis.axis_label = 'glucose uptake [mmol/(gDw.h)]'
        fig.yaxis.axis_label = ratio + ' [g/gDw]'

        if 'mrna' in ratio:
            fig.legend.location = 'bottom_right'

        output_file('plots/benchmark_{}.html'.format(ratio))
        show(fig)
        curdoc().clear()

        figures[ratio] = fig


    # mrna / ribosome content

    for ratio in ['mrna_ratio','prot_ratio']:
        fig = figure()
        for e, (model_name, df) in enumerate(model_data.items()):
            plot_prop_vs_growth(fig,df,ratio,model_name,cmap[e], plot_lines = False)

        fig.title.text = ratio + ' vs growth rate'
        fig.xaxis.axis_label = 'Growth rate [1/h]'
        fig.yaxis.axis_label = ratio + ' [g/gDw]'

        if 'mrna' in ratio:
            fig.legend.location = 'bottom_right'
            fig.square(neidhardt_mu, neidhardt_rrel,
                       line_color=cmap[3],
                       fill_alpha = 0,
                       size = 10,
                       line_width = 2,
                       legend='Neidhardt et al. measurements')
        elif 'prot' in ratio:
            fig.square(neidhardt_mu, neidhardt_prel,
                       line_color=cmap[3],
                       fill_alpha = 0,
                       size = 10,
                       line_width = 2,
                       legend='Neidhardt et al. measurements')

        output_file('plots/benchmark_{}_vs_mu.html'.format(ratio))
        show(fig)
        curdoc().clear()

        figures[ratio] = fig
