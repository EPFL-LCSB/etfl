from etfl.tests.small_model import create_etfl_model
from etfl.analysis.dynamic import run_dynamic_etfl
from etfl.optim.config import standard_solver_config
from etfl.core.reactions import EnzymaticReaction
from etfl.optim.constraints import ForwardCatalyticConstraint, \
    BackwardCatalyticConstraint, ExpressionCoupling
from etfl.io.json import load_json_model

from pytfa.utils.logger import get_timestr

import bokeh.plotting as bp
from bokeh.palettes import Category10
from bokeh.layouts import row, gridplot
from bokeh.models import LinearAxis, Range1d

from math import exp
import pandas as pd

from os.path import join
from os import makedirs

import yaml

def read_config(yaml_file):
    with open(yaml_file, 'rb') as f:
        conf = yaml.load(f.read(), Loader=yaml.SafeLoader)  # load the config file
    return conf


def run_detfl(yaml_file, uptake_funs, medium_funs=dict(), model=None, ini_sol=None):    # pass in

    params = read_config(yaml_file)

    for key, value in params.items():
        print('%s: %s' % (key, value))

    if model is None:
        model = load_json_model(params['model'])


    standard_solver_config(model)
    model.solver.configuration.verbosity = params['options']['verbose']

    time_data = run_dynamic_etfl(model,
                                 timestep=params['simulation']['timestep'],
                                 tfinal=params['simulation']['tfinal'],
                                 uptake_fun=uptake_funs,
                                 medium_fun=medium_funs,
                                 uptake_enz=params['assumptions']['uptake_enz'],
                                 S0=params['assumptions']['S0'],
                                 X0=params['assumptions']['X0'],
                                 inplace=params['options']['inplace'],
                                 initial_solution = ini_sol,
                                 chebyshev_include = params['options']['chebyshev_include'],
                                 dynamic_constraints = params['options']['constraints']
                                 )

    # plot_dynamics(model, time_data)
    out_path = join('outputs',params['tag']+get_timestr())
    makedirs(out_path)
    time_data.to_csv(join(out_path,'solution.csv'))
    write_yaml(params,join(out_path,'config.yaml'))

    return time_data

def write_yaml(config, path):
    with open(path, 'w') as outfile:
        yaml.dump(config, outfile)

def get_enzymes_of_subsystem(model, subsystem):
    reactions = [x for x in model.reactions if subsystem.lower() in x.subsystem.lower()]

    enzymes = [item for rxn in reactions if hasattr(rxn,'enzymes') and rxn.enzymes is not None
               for item in rxn.enzymes
               ]

    return enzymes
