from plotting import summarize_model
import pandas as pd
import sys
from os.path import join
from detfl_wrapper import read_config


import argparse

parser = argparse.ArgumentParser(description='Plots dETFL results')
parser.add_argument('output_path', nargs='?',
                    help='path to solution.csv and config.yaml')
parser.add_argument('--backend', type=str, nargs='?', default='canvas',
                    help='canvas (png) or svg')
args = parser.parse_args()

# model = None

fluxes = ['EX_glc__D_e','EX_lcts_e','EX_o2_e','EX_ac_e']

glc_fluxes  = ['GLCDpp','GLCabcpp','GLCptspp','GLCt2pp','total']
lcts_fluxes = ['LACZpp','LCTSt3ipp','LCTStpp','total']
ac_fluxes   = ['ACt2rpp', 'ACt4pp', 'total']

# lcts_enzymes = ['EZ_'+ x for x in ['LACZpp_EG12013_MONOMER',
#                                    'LCTSt3ipp_B2170_MONOMER',
#                                    'LCTSt3ipp_B0070_MONOMER',
#                                    'LCTSt3ipp_YDEA_MONOMER',
#                                    'LCTStpp_LACY_MONOMER',
#                                   ]] + ['total']
# glc_enzymes = ['EZ_'+ x for x in [
#                                   'GLCDpp_GLUCDEHYDROG_MONOMER_mod_pqq',
#                                   'GLCDpp_G6437_MONOMER_mod_ca2_mod_pqq',
#                                   'GLCabcpp_ABC_18_CPLX',
#                                   'GLCptspp_CPLX_164',
#                                   'GLCptspp_CPLX_157',
#                                   'GLCt2pp_GALP_MONOMER',
#                                   ]] + ['total']
# ac_enzymes = ['EZ_' + x for x in ['ACt4pp_YJCG_DASH_MONOMER',]]

ac_enzymes   = ['EZ_'+ x for x in ['ACt4pp_b4067',
                                   'ACt2rpp_b0010',
                                  ]] + ['total']
lcts_enzymes = ['EZ_'+ x for x in [
                                   'LACZpp_EG12013_MONOMER',
                                   'LCTSt3ipp_B2170_MONOMER',
                                   'LCTSt3ipp_B0070_MONOMER',
                                   'LCTSt3ipp_YDEA_MONOMER',
                                   'LCTStpp_LACY_MONOMER',
                                  ]] + ['total']
glc_enzymes = ['EZ_'+ x for x in [
                                  'GLCDpp_GLUCDEHYDROG_MONOMER_mod_pqq',
                                  'GLCDpp_G6437_MONOMER_mod_ca2_mod_pqq',
                                  'GLCabcpp_ABC_18_CPLX',
                                  'GLCptspp_CPLX_164',
                                  'GLCptspp_CPLX_157',
                                  'GLCt2pp_GALP_MONOMER',
                                  ]] + ['total']


species = ['S_EX_' + x for x in ['glc__D_e','lcts_e','o2_e','ac_e']]
groups = {'fluxes':fluxes,
          'glc_enzymes':glc_enzymes,
          'lcts_enzymes':lcts_enzymes,
          'ac_enzymes':ac_enzymes,
          'glc_fluxes':glc_fluxes,
          'lcts_fluxes':lcts_fluxes,
          'ac_fluxes':ac_fluxes,
          'species':species}

if __name__ == '__main__':

    if not args.output_path:
        time_data_path = 'tmp_detfl.csv'
        model_tag = 'tmp'
        output_path = 'plots'
    else:
        output_path = args.output_path
        config = read_config(join(output_path,'config.yaml'))
        model_tag = config['tag']
        time_data_path = join(output_path,'solution.csv')

    time_data = pd.read_csv(time_data_path, header=0, index_col=0)

    summarize_model(None, time_data, groups, output_path=output_path,
                    model_tag=model_tag, backend=args.backend)
