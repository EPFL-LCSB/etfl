import bokeh.plotting as bp
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, LinearAxis, Range1d
from bokeh.palettes import viridis
from bokeh.transform import linear_cmap
import pandas as pd

from os.path import join

data_folder = 'data'
plots_folder = 'plots'
my_data_folder = join('outputs','tmp')

T_OFFSET_ENJ2VAR = 4
T_OFFSET_ME2VAR = 2

def get_varma_data():
    varma_ac  = pd.read_csv(join(data_folder, 'varma7_ac.csv'), header=0)
    varma_ac. columns = ['t', 'ac']

    varma_glc = pd.read_csv(join(data_folder, 'varma7_gl.csv'), header=0)
    varma_glc.columns = ['t', 'glc']

    varma_bio = pd.read_csv(join(data_folder, 'varma7_bm.csv'), header=0)
    varma_bio.columns = ['t', 'biomass']

    centered_t = pd.concat([varma_bio['t'],varma_glc['t'],varma_ac['t']],axis=1).mean(axis=1)

    centered_dat = pd.DataFrame(columns = ['t','ac','glc','X'])

    centered_dat['t']   = centered_t
    centered_dat['ac']  = varma_ac['ac']
    centered_dat['glc'] = varma_glc['glc']
    centered_dat['biomass'] = varma_bio['biomass']

    return centered_dat

def get_enjalbert_data():
    raw_dat  = pd.read_csv(join(data_folder, 'fig2_fromEnjalbert2015.csv'), header=0)
    raw_dat.columns = ['tmin','t','od600','od600_std','glc','glc_std','ac','ac_std']
    raw_dat['t'] += T_OFFSET_ENJ2VAR
    return raw_dat

def get_my_data():
    # raw_data = pd.read_csv(join(my_data_folder,'solution.csv'),header=0,index_col=0)
    raw_data = pd.read_csv(join(my_data_folder,'tmp_detfl.csv'),header=0,index_col=0)
    out = raw_data.loc[['t','S_EX_glc__D_e','S_EX_ac_e','X']].T
    out['t'] += T_OFFSET_ME2VAR
    return out

varma = get_varma_data()
enjalbert = get_enjalbert_data()
detfl = get_my_data()

bp.output_file(join(plots_folder,'exp_acetate.html'))

def plot_varma(p,x,y,**kwargs):
    p.circle(x=x,y=y,size=7,line_width=3,fill_alpha=0, **kwargs)

def plot_enjalbert(p,x,y,**kwargs):
    p.square(x=x,y=y,size=7,line_width=3,fill_alpha=0,color='green',**kwargs)

def plot_me(p,x,y,**kwargs):
    p.cross(x=x,y=y,size=7,line_width=3,fill_alpha=0,color='red',**kwargs)

p1 = bp.figure()

plot_varma    (p1, x=varma['t']    ,y=varma['ac']    )
plot_enjalbert(p1, x=enjalbert['t'],y=enjalbert['ac'])
plot_me       (p1, x=detfl['t']    ,y=detfl['S_EX_ac_e'])
p1.xaxis.axis_label = 't [h]'
p1.yaxis.axis_label = 'Acetate [mmol/(gDW)]'

p2 = bp.figure()

plot_varma    (p2, x=varma['t']    ,y=varma['glc']    )
plot_enjalbert(p2, x=enjalbert['t'],y=enjalbert['glc'])
plot_me       (p2, x=detfl['t']    ,y=detfl['S_EX_glc__D_e'])
p2.xaxis.axis_label = 't [h]'
p2.yaxis.axis_label = 'Glucose [mmol/(gDW)]'

p3 = bp.figure()

plot_varma(p3,x=varma['t'],y=varma['biomass'])
plot_me   (p3,x=detfl['t'],y=detfl['X'])

p3.y_range.start = 0
p3.y_range.end = 1.05 * varma['biomass'].max()

# Setting the second y axis range name and range
p3.extra_y_ranges = {"enjalbert": Range1d(start=0, end=1.05 * enjalbert['od600'].max())}

# Adding the second axis to the plot.
p3.add_layout(LinearAxis(y_range_name="enjalbert",
                       axis_label='OD600'), 'right')

plot_enjalbert(p3, x=enjalbert['t'],y=enjalbert['od600'],y_range_name='enjalbert')
p3.xaxis.axis_label = 't [h]'
p3.yaxis.axis_label = 'Cells [g/L]'

bp.show(column([p1,p2,p3]))