import os
import pandas as pd
import re
import numpy as np
import bokeh.plotting as bp

log_folder = './logs/'

OPTIM_REGEX = re.compile(r'optimize')
TIME_REGEX = re.compile(r'(\d*\.\d*) sec')

bp.curdoc().clear()

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

    measured = [item for subl in times.values() for item in subl]
    bins = np.logspace(np.log10(min(measured)),np.log10(max(measured)), 50)

    hist, edges = np.histogram(measured, density=False, bins=bins)

    p1 = bp.figure(title='Distribution of optimization times', x_axis_type='log')

    p1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
            fill_color="#036564", line_color="#033649")

    p1.legend.location = "center_right"
    p1.legend.background_fill_color = "darkgrey"
    p1.xaxis.axis_label = 'seconds'
    p1.yaxis.axis_label = 'counts'

    bp.output_file('./plots/timing_hist.html')
    bp.show(p1)



