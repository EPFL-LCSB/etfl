import os
import pandas as pd
import re
import numpy as np
from scipy.stats.mstats import gmean
import bokeh.plotting as bp
from collections import defaultdict
from bokeh.palettes import Category10
from bokeh.layouts import gridplot, column
from bokeh.models import ColumnDataSource, FixedTicker, PrintfTickFormatter, \
    Label

from scipy.stats.kde import gaussian_kde

cmap = Category10[10]

log_folder = './logs/'

OPTIM_REGEX = re.compile(r'optimize')
TIME_REGEX = re.compile(r'(\d*\.\d*) sec')
MODELTYPE_REGEX=re.compile(r'T([01])E[01]N([01])')

bool2name={
    (False, False)  : 'T0E1N0',
    (False, True)   : 'T0E1N1',
    (True, False)   : 'T1E1N0',
    (True, True)    : 'T1E1N1',
}

bp.curdoc().clear()

str2bool = lambda x:bool(int(x))

def get_file_optim_times(filename):

    times = []

    with open(os.path.join(log_folder,filename), 'r') as fid:
        for line in fid:
            if OPTIM_REGEX.search(line):
                times.append(float(TIME_REGEX.search(line).group(1)))

    return times

def logx_kde(x, y):
    values = np.log10(y)
    f_hat = gaussian_kde(values, np.std(values)/3)
    return np.power(10, f_hat(x))

def plot_hist(measured):
    all_times = [item for subl in measured.values() for item in subl]
    bins = np.logspace(np.log10(min(all_times)), np.log10(max(all_times)), 50)

    p = list()

    for e, (model_type, these_times) in enumerate(measured.items()):

        avg = np.mean(these_times)
        geo = gmean(these_times)
        med = np.median(these_times)
        pc05 = np.percentile(these_times, 5)
        pc95 = np.percentile(these_times,95)

        p1 = bp.figure(title='Distribution of optimization times \n'
                             '({} counts, average {:.2f} s)'
                       .format(len(these_times), np.mean(these_times)),
                       x_axis_type='log')
        hist, edges = np.histogram(these_times,
                                   density=False,
                                   # density=True,
                                   bins=bins)

        df_data = pd.DataFrame({'h':hist,
                                'left':edges[:-1],
                                'right':edges[1:],
                                })
        df_data['color'] = cmap[e]
        df_data['alpha'] = 1
        df_data['alpha'].loc[df_data['right'] < pc05] = 0.5
        df_data['alpha'].loc[df_data['left'] > pc95] = 0.5

        cds = ColumnDataSource(df_data)

        p1.quad(top='h', bottom=0,
                left='left', right='right',
                fill_color='color',
                fill_alpha='alpha',
                line_color=None,
                source = cds,
                legend=model_type)
        ymax = max(hist)


        p1.line([avg,avg],[0, ymax], legend='Mean: {:.2f}s'.format(avg), line_color = 'black')
        p1.line([geo,geo],[0, ymax], line_dash='dashed', line_color = 'black',
                legend='Geometric Mean: {:.2f}s'.format(geo))
        p1.line([med,med],[0, ymax], line_dash='dotted', line_color = 'black',
                legend='Median: {:.2f}s'.format(med))

        # p1.add_layout(Label(x=avg,y=ymax,text='{:.2f}s'.format(avg)))
        # p1.add_layout(Label(x=med,y=ymax,text='{:.2f}s'.format(med)))
        # p1.add_layout(Label(x=geo,y=ymax,text='{:.2f}s'.format(geo)))

        # last_height += hist
        p1.legend.location = "center_right"
        p1.xaxis.axis_label = 'seconds'
        p1.yaxis.axis_label = 'counts'

        p.append(p1)

    return gridplot([[p[0],p[1]],[p[2],p[3]]])
    # return column(p)

def plot_ridge(measured):
    all_times = [item for subl in measured.values() for item in subl]

    x = np.logspace(np.log10(min(all_times))-1, np.log10(max(all_times))+1, 500)
    source = ColumnDataSource(data=dict(x=x))

    p = bp.figure(title='Distribution of optimization times ({} measurements)'
                   .format(len(all_times)),
                   x_axis_type='log')
    for i, (model_type, these_times) in enumerate(measured.items()):
        if len(these_times) == 0:
            continue
        y = logx_kde(x,these_times) + 1
        source.add(y, model_type)
        p.patch('x', model_type, color=cmap[i], alpha=0.6, line_color="black", source=source)


    p.y_range.range_padding = 0.12

    return p


if __name__ == '__main__':

    times = dict()

    for filename in os.listdir(log_folder):
        if filename.endswith(".log"):

            times[filename] = get_file_optim_times(filename)

        else:
            continue

    measured = defaultdict(list)

    for filename,these_t in times.items():
        model_type = MODELTYPE_REGEX.search(filename)
        if model_type:
            has_thermo, has_alloc = model_type.groups(1)
            measured[bool2name[str2bool(has_thermo),str2bool(has_alloc)]] += these_t
        else:
            pass
            # measured['other'] += these_t

    all_times = [item for subl in times.values() for item in subl]

    p1 = plot_hist(measured)
    # p2 = plot_ridge(measured)

    bp.output_file('./plots/timing_hist.html')
    bp.show(p1)
    # bp.show(p2)



