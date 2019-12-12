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
from etfl.core.equilibrium import add_rnap_binding_equilibrium_constraints
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

def write_conf(model, conf):
    filename = 'outputs/' + conf['tag'] + '/config.yaml'
    conf['objective'] = str(model.objective.expression)
    with open(filename,'w') as fid:
        yaml.dump(config, fid)



if __name__ == '__main__':
    config = read_config(CONFIG)

    if config['model'] != 'debug':
        model = load_json_model(config['model'],solver=config['options']['solver'])
    else:
        from etfl.tests.small_model import create_etfl_model
        model = create_etfl_model(0,1, solver=config['options']['solver'])

    standard_solver_config(model)
    model.solver.configuration.verbosity = config['options']['verbose']

    copy_number = int(config['simulation']['copy_number'])
    vector_generator = vec_dict[config['simulation']['vector']]
    has_rnap = config['simulation']['add_rnap']

    my_plasmid = vector_generator(model, has_rnap)


    #####################
    # Model integration #
    #####################

    # 1. Adding RNAP eq constraints
    # Constant for binding of RNAP to lac promoter
    # Value 	550 nM
    # Organism 	Bacteria Escherichia coli
    # Reference 	Bintu L, Buchler NE, Garcia HG, Gerland U, Hwa T, Kondev J, Phillips R. Transcriptional regulation by the numbers: models. Curr Opin Genet Dev. 2005 Apr15(2):116-24. Fig. 1PubMed ID15797194
    # Primary Source 	[38] Liu M, Gupte G, Roy S, Bandwar RP, Patel SS, Garges S. Kinetics of transcription initiation at lacP1. Multiple roles of cyclic AMP receptor protein. J Biol Chem. 2003 Oct 10 278(41):39755-61PubMed ID12881519
    # Method 	Primary source abstract: "[Researchers] examined the kinetics of open complex formation at the lacP1 promoter using tryptophan fluorescence of RNA polymerase and DNA fragments with 2-aminopurine substituted at specific positions."
    # Comments 	P.118 caption to fig.1c: "In particular, making the simplest assumption that the genomic background for RNAP is given only by the non-specific binding of RNAP with DNA, [investigators] take K[NS]pd=10 000 nM [ref 37], for the lac promoter K[S]pd=550 nM [primary source] and for the T7 promoter, K[S]pd=3nM [BNID 103592]. For the lac promoter, this results in Δɛpd=-2.9kBT and for the T7 promoter, Δɛpd=-8.1kBT."
    # BNID 	103590
    #          ecoli density is ~ 1.2
    #    mol/L * L/kg   * kg/g * g/gDW
    Kb = 550e-9 * 1/1.2 * 1000 * 1/0.5
    add_rnap_binding_equilibrium_constraints(model,
                                             the_rnap=model.rnap['rnap'],
                                             Kb=Kb)

    # Sigma factor
    from etfl.data.ecoli import get_sigma_70
    sigma, holo = get_sigma_70(model.rnap['rnap'])
    model.add_sigma_factor('rnap', sigma, holo)

    1/0

    transmodel = TransModel(model, inplace = config['options']['inplace'])
    transmodel.add_vector(my_plasmid, copy_number = copy_number)
    transmodel.reactions.ALS.upper_bound = 0

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
        #TODO: fix
    # write_conf(transmodel, config)
