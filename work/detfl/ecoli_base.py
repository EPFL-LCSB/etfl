from detfl_wrapper import run_detfl, read_config
from etfl.io.json import load_json_model
from etfl.optim.config import standard_solver_config
import yaml
from math import exp
import sys

if len(sys.argv)<=1:
    CONFIG = 'monod_acetate.yml'
    # CONFIG = 'monod_lactose.yml'
else:
    CONFIG = sys.argv[1]

print('Using configuration file: {}'.format(CONFIG))

def get_uptake_funs():

    uptake_funs = dict()

    # Glucose:

    Vmax0 = 10 # mmol/(h.mmol[E]) Mahadevan et al. 2002
    # Vmax0 = 15 # mmol/(h.mmol[E])
    # Vmax0 = 1
    Km0 = 0.015 # mM, Mahadevan et al. 2002, Wong et al. 1997

    uptake_funs['EX_glc__D_e'] = lambda x: Vmax0 * x / (Km0 + x)

    # Lactose:
    # Olsen, S. G., and R. J. Brooker.
    # "Analysis of the structural specificity of the lactose permease toward sugars."
    # Journal of Biological Chemistry 264.27 (1989): 15982-15987.
    # http://www.jbc.org/content/264/27/15982.full.pdf+html
    # TODO: Stop assuming 1mgProt/gDW
    # Vmax_lac = 210 /1000 * 60 *1  # nmol/min/mgProt * mmol/nmol * min/h * mgProt/gDW
    Vmax_lac = 1
    Km_lac = 1.3 # mM

    uptake_funs['EX_lcts_e'] = lambda x: Vmax_lac * x / (Km_lac + x)

    # O2
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC134846/
    # Alexeeva S, Hellingwerf KJ, Teixeira de Mattos MJ.
    # Quantitative assessment of oxygen availability: perceived aerobiosis and
    # its effect on flux distribution in the respiratory chain of
    # Escherichia coli.
    # J Bacteriol. 2002;184(5):1402-6.
    # The cytochrome bo oxidase has a Km for O2 of 2 × 10−4 mM and a Vmax of
    # 6.6 × 10−2 mmol of O2·nmol of cytochrome o−1·h−1 (18). This corresponds,
    # at the measured rDOT value of 1.6 × 10−2 mM, to a cytochrome bo oxidase
    # content of 73 nmol of protein·g (dry weight)−1
    # VmaxO2 = 6.6*1e-2 *73 # mmol/h
    VmaxO2 = 15 # mmol/h
    KmO2 = 2 * 1e-4 #mM
    uptake_funs['EX_o2_e'] = lambda x: VmaxO2 * x /(KmO2 + x)

    # Acetate
    fix_ac = lambda x: 15 if x>0 else 0
    succurro = lambda x: 10 * x /(0.01 + x)
    succurro_v2 = lambda x: 3 * x /(0.01 + x)
    uptake_funs['EX_ac_e'] = succurro_v2

    return uptake_funs

def prepare_model(in_model, config, uptake_fun):
    """
    :param in_model:
    :return:
    """

    #  ('LCTSt3ipp',
    # 'Lactose transport via proton aniport _periplasm',
    # 'h_p + lcts_c --> h_c + lcts_p'),
    # We set it to 0, since we will not produce lactose,
    # and do not have any enzyme related to it

    v0  = config['assumptions']['v0']
    S0  = config['assumptions']['S0']
    ubs = config['assumptions']['upper_bounds']

    try:
        model.reactions.LCTSt3ipp.lower_bound = 0
    except AttributeError:
        pass

    for rxn_id,lb in v0.items():
        try:
            model.reactions.get_by_id(rxn_id).lower_bound = lb
        except KeyError: # For debug models
            model.logger.warning('Reaction {} not in model - could not '
                                 'initialize flux'.format(rxn_id))

    for rxn_id,ub in ubs.items():
        try:
            model.reactions.get_by_id(rxn_id).upper_bound = ub
        except KeyError: # For debug models
            model.logger.warning('Reaction {} not in model - could not '
                                 'set upper bound'.format(rxn_id))


    sol_ini = in_model.optimize()

    return sol_ini


def get_medium_funs(config):


    try:
        mode = config['simulation']['mode']
    except KeyError:
        mode = 'batch'

    timestep = config['simulation']['timestep']

    epsilon = timestep / 100
    try:
        S0_o2 = config['assumptions']['S0']['EX_o2_e']  # mmol/L
        S1_o2 = 0.21  # mmol/L
    except KeyError:
        print('No o2 S0 found')
        S0_o2 = 0
        S1_o2 = 0
    kla_o2 = 7.5  # h^-1

    S0_glc = config['assumptions']['S0']['EX_glc__D_e']  # mmol/L
    S1_glc = 1  # mmol/L

    try:
        S0_lac = config['assumptions']['S0']['EX_lcts_e']  # mmol/L
        S1_lac = 10  # mmol/L
    except KeyError:
        print('No lactose S0 found')
        S0_lac = 0
        S1_lac = 0

    S0_ac = config['assumptions']['S0']['EX_ac_e']

    X0 = config['assumptions']['X0']

    glc_fun = lambda t, S, S0=S0_glc, S1=S1_glc: \
        max(S, 0)

    glc_constant = lambda t, S, S0=S0_glc: S0
    glc_switch1  = lambda t, S, S0=S1_glc: S0 if t > 0.5 else 0

    # S + S1 if abs(t - 1) <= timestep+epsilon and S <= S0 else max(S, 0)
    # S1 if t > 1  else S0

    lac_fun = lambda t, S, S0=S0_lac, S1=S1_lac: \
        max(S, 0)

    lac_constant = lambda t, S, S0=S0_lac, S1=S1_lac: S0
    lac_switch1  = lambda t, S, S0=S0_lac: max(S,0) if t > 0.5 else S0


    o2_fun = lambda t, S, S0=S0_o2, S1=S1_o2: \
        S0  # if t > 1  else S0
    # S + S1 if abs(t - 1) <= timestep+epsilon and S <= S0 else max(S, 0)

    glc_free = lambda t, S: max(S0_glc, 0)

    # Integrated linearization of the diffusion over dt
    # o2_diff = lambda t, S, S0=S0_o2, kla=kla_o2: max(S0 - (S0 - S) * exp(-kla * timestep), 0)
    o2_diff = lambda t, S, S0=S0_o2, kla=kla_o2: max(S + kla * (S0 - S) * timestep, 0)
    o2_constant = lambda t, S, S0=S0_o2: S0

    ac_fun = lambda t, S: max(S, 0)

    if mode == 'batch':
        medium_funs = {
            'EX_glc__D_e': glc_fun,
            'EX_lcts_e': lac_fun,
            'EX_o2_e': o2_diff,
            'EX_ac_e': ac_fun,
        }
    elif mode == 'chemostat_lcts':
        medium_funs = {
            'EX_glc__D_e': glc_fun,
            'EX_lcts_e': lac_constant,
            'EX_o2_e': o2_constant,
            'EX_ac_e': ac_fun,
        }
    elif mode == 'chemostat_glc':
        medium_funs = {
            'EX_glc__D_e': glc_constant,
            'EX_lcts_e': lac_fun,
            'EX_o2_e': o2_constant,
            'EX_ac_e': ac_fun,
        }
    elif mode == 'switch1':
        medium_funs = {
            'EX_glc__D_e': glc_switch1,
            'EX_lcts_e': lac_switch1,
            'EX_o2_e': o2_constant,
            'EX_ac_e': ac_fun,
        }
    elif mode == 'succurro':
        medium_funs = {
            'EX_glc__D_e': glc_fun,
            'EX_lcts_e': lac_fun,
            'EX_o2_e': o2_constant,
            'EX_ac_e': ac_fun,
        }
    elif mode == 'varma':
        medium_funs = {
            'EX_glc__D_e': glc_fun,
            'EX_lcts_e': lac_fun,
            'EX_o2_e': o2_constant,
            'EX_ac_e': ac_fun,
        }


    return  medium_funs


if __name__ == '__main__':
    config = read_config(CONFIG)

    has_lcts = 'EX_lcts_e' in config['assumptions']['S0']
    has_o2   = 'EX_o2_e'   in config['assumptions']['S0']

    medium_funs = get_medium_funs(config)
    uptake_funs = get_uptake_funs()

    if not has_lcts:
        medium_funs.pop('EX_lcts_e')
        uptake_funs.pop('EX_lcts_e')

    if not has_o2:
        medium_funs.pop('EX_o2_e')
        uptake_funs.pop('EX_o2_e')

    if config['model'] != 'debug':
        model = load_json_model(config['model'])
    else:
        from etfl.tests.small_model import create_etfl_model
        model = create_etfl_model(0,0)

    standard_solver_config(model)
    model.solver.configuration.verbosity = 0

    ini_sol = prepare_model(model,
                            config,
                            uptake_fun = uptake_funs)

    time_data = run_detfl(model=model,
                          yaml_file=CONFIG,
                          ini_sol=ini_sol,
                          uptake_funs=uptake_funs,
                          medium_funs=medium_funs)

