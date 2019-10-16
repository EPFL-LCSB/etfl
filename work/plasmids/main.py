# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Models vectors in ETFL models

.. moduleauthor:: ETFL team

Model RNAP limitation in case of vector addition to a host
"""

from etfl.io.json import load_json_model
from etfl.optim.config import standard_solver_config
from etfl.analysis.summary import print_standard_sol
from etfl.core.vector import TransModel
import yaml
import sys

from utils import read_config
from vector import get_bdo_plasmid, get_debug_plasmid
import analysis

if len(sys.argv)<=1:
    raise ValueError('Please provide a YAML configuration file')
else:
    CONFIG = sys.argv[1]

print('Using configuration file: {}'.format(CONFIG))

vec_dict={'plasmid_pET-AR-ALS':get_bdo_plasmid,
          'debug':get_debug_plasmid}

if __name__ == '__main__':
    config = read_config(CONFIG)

    if config['model'] != 'debug':
        model = load_json_model(config['model'])
    else:
        from etfl.tests.small_model import create_etfl_model
        model = create_etfl_model(0,1)

    standard_solver_config(model)
    model.solver.configuration.verbosity = config['options']['verbose']

    copy_number = config['simulation']['copy_number']
    vector_generator = vec_dict[config['simulation']['vector']]
    has_rnap = config['simulation']['add_rnap']

    my_plasmid = vector_generator(model, has_rnap)


    #####################
    # Model integration #
    #####################

    transmodel = TransModel(model, inplace = config['options']['inplace'])
    transmodel.add_vector(my_plasmid, copy_number = copy_number)

    transmodel.optimize()
    print_standard_sol(transmodel)

    outputs = dict()

    if 'chebyshev' in config['analysis']:
        arguments = config['analysis'].pop('chebyshev')
        arguments['model'] = transmodel
        # arguments['tag'] = config['tag']
        outputs['chebyshev'] = analysis.chebyshev(**arguments)

    for analysis_str,arguments in config['analysis'].items():
        if arguments != False:
            if arguments == True:
                arguments = dict()
            arguments['model'] = transmodel
            arguments['tag'] = config['tag']
            outputs[analysis_str] = getattr(analysis, analysis_str)(**arguments)

        # Check antibiotic resistance is expressed ?
