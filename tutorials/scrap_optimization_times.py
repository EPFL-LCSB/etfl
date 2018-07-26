import os
import pandas as pd
import re
import numpy as np
import bokeh.plotting as bp
from collections import defaultdict
from bokeh.palettes import Category10

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
            measured['other'] += these_t

    all_times = [item for subl in times.values() for item in subl]
    bins = np.logspace(np.log10(min(all_times)), np.log10(max(all_times)), 50)


    p1 = bp.figure(title='Distribution of optimization times ({} measurements)'
                   .format(len(all_times)),
                   x_axis_type='log')

    last_height = np.array([0]*(len(bins)-1))

    for e, (model_type, these_times) in enumerate(measured.items()):

        hist, edges = np.histogram(these_times, density=False, bins=bins)

        p1.quad(top=hist+last_height, bottom=last_height,
                left=edges[:-1], right=edges[1:],
                fill_color=cmap[e],
                # fill_alpha=0.75,
                line_color=None,
                legend=model_type)

        last_height += hist

    p1.legend.location = "center_right"
    p1.legend.background_fill_color = "darkgrey"
    p1.xaxis.axis_label = 'seconds'
    p1.yaxis.axis_label = 'counts'

    bp.output_file('./plots/timing_hist.html')
    bp.show(p1)



