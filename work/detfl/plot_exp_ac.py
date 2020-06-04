import bokeh.plotting as bp
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, LinearAxis, Range1d
from bokeh.palettes import viridis
from bokeh.transform import linear_cmap
from bokeh.io import export_svgs
import pandas as pd
import numpy as np

from os.path import join
from os import makedirs

from scipy import signal

data_folder = 'data'
plots_folder = 'plots'
subplots_folder = join(plots_folder, 'exp_acetate_sub')
makedirs(subplots_folder,exist_ok=True)

# my_data_folder_enj   = join('outputs', 'vETFL_enjalbert_vmax_20200330_134009') #1
# my_data_folder_enj   = join('outputs', 'vETFL_enjalbert_vmax_20200330_083510') #2
# my_data_folder_enj   = join('outputs', 'vETFL_enjalbert_20200407_210513') #new
my_data_folder_enj   = join('outputs', 'vETFL_enjalbert_20200414_182924') #new
# my_data_folder_varma = join('outputs', 'vETFL_varma_20200331_224535') #1
# my_data_folder_varma = join('outputs', 'vETFL_varma_vmax_20200330_131840') #2
# my_data_folder_varma = join('outputs', 'vETFL_varma_20200407_190236') #new
my_data_folder_varma = join('outputs', 'vETFL_varma_20200414_152012') #new

BACKEND = 'svg'

def _resample(t,x,delta_t,num=50):
    new_t = np.arange(t.min(),t.max(),delta_t)
    x_hat,t_hat = signal.resample(x,num,new_t)
    return x_hat,t_hat

def _get_offset(source,target,y_col,target_y_col='biomass'):
    #TODO: make time-aware
    dx = np.mean(np.diff(source.t.values))
    source_hat,_ = _resample(source.t,source[y_col],dx)
    target_hat,_ = _resample(target.t,target[target_y_col],dx)
    shift = (np.argmax(signal.correlate(source_hat, target_hat)) - len(target_hat)) * dx
    return shift

def get_avg_offset(source,target,biomass_col='biomass'):
    o_ac  = _get_offset(source,target,'ac')
    o_glc = _get_offset(source,target,'glc')
    o_bio = _get_offset(source,target,biomass_col)

    return np.mean([o_ac,o_glc,o_bio])


def get_varma_data():
    varma_ac  = pd.read_csv(join(data_folder, 'varma7_ac.csv'), header=0)
    varma_ac. columns = ['t', 'ac']
    varma_ac = varma_ac.iloc[varma_ac['t'].sort_values().index]

    varma_glc = pd.read_csv(join(data_folder, 'varma7_gl.csv'), header=0)
    varma_glc.columns = ['t', 'glc']
    varma_glc = varma_glc.iloc[varma_glc['t'].sort_values().index]


    varma_bio = pd.read_csv(join(data_folder, 'varma7_bm.csv'), header=0)
    varma_bio.columns = ['t', 'biomass']
    varma_bio = varma_bio.iloc[varma_bio['t'].sort_values().index]


    centered_t = pd.concat([varma_bio['t'],varma_glc['t'],varma_ac['t']],axis=1).mean(axis=1)

    centered_dat = pd.DataFrame(columns = ['t','ac','glc','biomass'])

    centered_dat['t']   = centered_t
    centered_dat['ac']  = varma_ac['ac']
    centered_dat['glc'] = varma_glc['glc']
    centered_dat['biomass'] = varma_bio['biomass']

    centered_dat['t'] -= centered_dat[centered_dat['glc']<=centered_dat['glc'].max()*0.37]['t'].iloc[0]

    return centered_dat

def get_enjalbert_data():
    raw_dat  = pd.read_csv(join(data_folder, 'fig2_fromEnjalbert2015.csv'), header=0)
    raw_dat.columns = ['tmin','t','od600','od600_std','glc','glc_std','ac','ac_std']
    # raw_dat['t'] -= get_avg_offset(raw_dat,varma,'od600')
    raw_dat['t'] -= raw_dat[raw_dat['glc']<=raw_dat['glc'].max()*0.37]['t'].iloc[0]
    return raw_dat

def get_my_varma():
    raw_data = pd.read_csv(join(my_data_folder_varma, 'solution.csv'), header=0, index_col=0)
    # raw_data = pd.read_csv(join('tmp_detfl.csv'),header=0,index_col=0)
    # raw_data = pd.read_csv(join('outputs','tmp','tmp_detfl.csv'),header=0,index_col=0)
    out = raw_data.loc[['t','S_EX_glc__D_e','S_EX_ac_e','X','mu']].T
    out.columns = ['t','glc','ac','X','mu']
    # out['t'] -= get_avg_offset(out,varma,'X')
    out['t'] -= out[out['glc']<=out['glc'].max()*0.37]['t'].iloc[0]
    return out

def get_my_enj():
    raw_data = pd.read_csv(join(my_data_folder_enj, 'solution.csv'), header=0, index_col=0)
    # raw_data = pd.read_csv(join('outputs','tmp','tmp_detfl.csv'),header=0,index_col=0)
    out = raw_data.loc[['t','S_EX_glc__D_e','S_EX_ac_e','X','mu']].T
    out.columns = ['t','glc','ac','X','mu']
    # out['t'] -= get_avg_offset(out,varma,'X')
    out['t'] -= out[out['glc']<=out['glc'].max()*0.37]['t'].iloc[0]
    return out

varma     = get_varma_data()
enjalbert = get_enjalbert_data()
my_enj    = get_my_enj()
my_varma  = get_my_varma()

bp.output_file(join(plots_folder,'exp_acetate.html'))

def plot_varma(p,x,y,**kwargs):
    p.circle(x=x,y=y,size=7,line_width=3,fill_alpha=0, **kwargs)

def plot_enjalbert(p,x,y,**kwargs):
    p.square(x=x,y=y,size=7,line_width=3,fill_alpha=0,color='green',**kwargs)

def plot_my_enj(p, x, y, **kwargs):
    p.cross(x=x,y=y,size=7,line_width=3,fill_alpha=0,color='red',**kwargs)

def plot_my_varma(p, x, y, **kwargs):
    p.x(x=x,y=y,size=7,line_width=3,fill_alpha=0,color='red',**kwargs)

def plot_exp(p, x, y, **kwargs):
    p.x(x=x, y=y, size=7, line_width=2, fill_alpha=0, color='black', **kwargs)
def plot_sim(p,x,y,**kwargs):
    p.line(x=x,y=y,line_width=3, color='black',**kwargs)

def plot_conc(exp_dataset, sim_dataset, exp_label):

    p_ac = bp.figure(height = 300, width = 800)
    plot_exp(p_ac, x=exp_dataset['t'], y=exp_dataset['ac'],legend_label=exp_label + ' dataset')
    plot_sim(p_ac, x=sim_dataset['t'], y=sim_dataset['ac'],legend_label='simulated')
    p_ac.xaxis.axis_label = 't [h]'
    p_ac.yaxis.axis_label = 'Acetate [mmol/(L)]'

    p_glc = bp.figure(height = 300, width = 800)
    plot_exp(p_glc, x=exp_dataset['t'], y=exp_dataset['glc'],legend_label=exp_label + ' dataset')
    plot_sim(p_glc, x=sim_dataset['t'], y=sim_dataset['glc'],legend_label='simulated')
    p_glc.xaxis.axis_label = 't [h]'
    p_glc.yaxis.axis_label = 'Glucose [mmol/(L)]'

    p_ac .output_backend = BACKEND
    p_glc.output_backend = BACKEND

    return p_ac,p_glc

def plot_biomass(exp_dataset, sim_dataset, exp_label, biomass_col, is_OD = False):
    from plotting import make_growth_plot

    p = make_growth_plot(None,sim_dataset.T)

    if is_OD:
        # Setting the second y axis range name and range
        # p.extra_y_ranges = {'od_range': Range1d(start=0, end=1.05 * exp_dataset[biomass_col].max())}
        # # Adding the second axis to the plot.
        # p.add_layout(LinearAxis(y_range_name='od_range',
        #                          axis_label='OD'), 'right')
        # plot_exp(p, x=exp_dataset['t'], y=exp_dataset[biomass_col],legend_label=exp_label+' dataset',y_range_name='od_range')
        # p.y_range.start = 0
        # p.y_range.end = 1.05 * exp_dataset[biomass_col].max()

        # Instead, we rescale the OD:
        new_y = exp_dataset[biomass_col]
        new_y *= (sim_dataset['X'].max() - sim_dataset['X'].min()) / (new_y.max() - new_y.min())
        plot_exp(p, x=exp_dataset['t'], y=new_y                   ,legend_label=exp_label + ' dataset (rescaled)')
    else:
        plot_exp(p, x=exp_dataset['t'], y=exp_dataset[biomass_col],legend_label=exp_label + ' dataset')

    p.output_backend = BACKEND

    return p


p_ac_varma, p_glc_varma = plot_conc(varma,my_varma,'Varma')
p_ac_enj, p_glc_enj = plot_conc(enjalbert,my_enj,'Enjalbert')

p_bio_varma = plot_biomass(varma    ,my_varma,'Varma'    ,biomass_col='biomass',is_OD=False)
p_bio_enj   = plot_biomass(enjalbert,my_enj  ,'Enjalbert',biomass_col='od600'  ,is_OD=True )

to_show = column([p_ac_varma, p_glc_varma, p_ac_enj, p_glc_enj, p_bio_varma, p_bio_enj])
export_svgs(to_show, filename = join(subplots_folder,'plot.svg'))

bp.show(to_show)