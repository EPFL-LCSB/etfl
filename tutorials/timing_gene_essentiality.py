import cobra.io.json
import pytfa.io.json
from etfl.io.json import load_json_model
from etfl.optim.config import standard_solver_config, gene_ko_config
import optlang

import numpy as np

from cobra.flux_analysis import single_gene_deletion

from tqdm import tqdm

import time

solver = 'optlang-gurobi'
USE_GPR = False
#
# ecoli_fba = cobra.io.json.load_json_model('models/iJO1366_T0E0N0__20180606_121758.json')
# ecoli_fba.solver = solver
# ecoli_tfa = pytfa.io.json.load_json_model('models/iJO1366_T1E0N0__20180606_121751.json')
# ecoli_tfa.solver = solver

# vETFL
# ecoli = load_json_model('models/SlackModel iJO1366_vETFL_v_0.12_431_enz_128_bins__20200605_083556.json')
# ecoli = load_json_model('models/SlackModel iJO1366_vETFL_mean_kcat_431_enz_128_bins__20200605_090017.json')
# vETFL_infer:
ecoli = load_json_model('models/SlackModel iJO1366_vETFL_infer_infer_2084_enz_128_bins__20200605_102233.json')
ecoli.solver = solver

standard_solver_config(ecoli)
ecoli.optimize()
max_growth = ecoli.solution.objective_value

print('Growth               : {}'.format(ecoli.solution.objective_value))
print(' - Ribosomes produced: {}'.format(ecoli.solution.raw.EZ_rib))
print(' - RNAP produced: {}'.format(ecoli.solution.raw.EZ_rnap))


# --


def safe_optim(model):
    try:
        ret = model.slim_optimize()
    except KeyError:
        ret = np.nan
    return ret


def ko_gene(model, gene_id):
    try:
        the_trans = model.get_translation(gene_id)
    except KeyError:
        return None

    try:
        the_trans.upper_bound = 0
        growth = model.slim_optimize()
        #
        # the_trans.upper_bound = initial_value
    except KeyError:
        growth = np.nan

    return growth


def timeit(method):
    def timed(*args, **kwargs):
        if 'log_data' in kwargs:
            log_data = kwargs.pop('log_data')
        else:
            log_data = dict()

        ti = time.time()
        result = method(*args, **kwargs)
        tf = time.time()


        name = method.__name__
        log_data[name] = tf-ti

        print('{}: {}s'.format(method.__name__, tf-ti))

        return result
    return timed

def gene_has_associated_enzyme(model, gene):
    return any([gene.id in x.composition for x in model.enzymes])

@timeit
def ko_etfl(model):
    growth = dict()
    gene_ko_config(model)
    model.solver.configuration.timeout = 14400

    model.growth_reaction.lower_bound = 0.1*max_growth
    # Gene essentiality is a feasibility problem
    # /!\ Also, integer 0 will get the reaction number 0 !
    # We need to use optlang's symbolic 0
    model.objective = optlang.symbolics.Zero

    for g in tqdm(model.genes):
        with model as model:

            if gene_has_associated_enzyme(model, g):
                ret = ko_gene(model, g.id)
                if ret is None:
                    # There is no translation associated to this gene
                    model.logger.warning('There is no translation associated '
                                         'to this gene: {}'.format(g.id))
                    continue
            else:
                # default to GPR KO
                if not USE_GPR:
                    continue # skip the gene
                model.genes.get_by_id(g.id).knock_out()
                ret = safe_optim(model)

            growth[g.id] = ret

    model.growth_reaction.lower_bound = 0
    return growth

@timeit
def ko_fba(model):
    return single_gene_deletion(model)

@timeit
def ko_tfa(model):
    return single_gene_deletion(model)


log_data = dict()

# ko_fba(ecoli_fba, log_data = log_data)
# ko_tfa(ecoli_tfa, log_data = log_data)
values = ko_etfl(ecoli, log_data = log_data)

import pandas as pd
df = pd.DataFrame.from_dict(values, orient='index')
filename = 'gene_essentiality_vETFL'
filename += '_gpr' if USE_GPR else ''
filename += '_mean' if '_mean' in ecoli.name else ''
filename += '_infer' if '_infer' in ecoli.name else ''

df.to_csv('outputs/{}.csv'.format(filename))
