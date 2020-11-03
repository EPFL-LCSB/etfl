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
gal_fluxes  = ['GALt2pp','GALabcpp','total']
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
glc_enzymes = ['EZ_'+ x for x in ['GLCDpp_GLUCDEHYDROG_MONOMER_mod_pqq',
                                  'GLCDpp_G6437_MONOMER_mod_ca2_mod_pqq',
                                  'GLCabcpp_ABC_18_CPLX',
                                  'GLCptspp_CPLX_157',
                                  'GLCptspp_CPLX_164',
                                  'GLCptspp_CPLX_165',
                                  'GLCt2pp_GALP_MONOMER',
                                  ]] + ['total']

gal_enzymes = ['EZ_' + x for x in [
                                  'GALt2pp_GALP_MONOMER',
                                  'GALabcpp_ABC_18_CPLX',
                                  'GALabcpp_ABC_46_CPLX',
                                  'GALKr_G7096_MONOMER',
                                  'GALKr_GALACTOKIN_MONOMER_mod_mg2',
                                  ]] + ['total']

species = ['S_EX_' + x for x in ['glc__D_e','lcts_e','o2_e','ac_e']]

# Checkpoint_fluxes
ackr_checkpoint = ['EZ_' + x for x in ['ACKr_ACETATEKINA_MONOMER_mod_mg2_1',
                                       'ACKr_GARTRANSFORMYL2_MONOMER_2',
                                       'ACKr_GARTRANSFORMYL2_MONOMER_2']] + ['total']

pyk_checkpoint = ['EZ_' + x for x in ['PYK_PKI_COMPLEX_mod_mg2_mod_k_2',
                                       'PYK_PKI_COMPLEX_mod_mn2_mod_k_3',
                                       'PYK_PKII_CPLX_mod_mg2_mod_k_0',
                                       'PYK_PKII_CPLX_mod_mg2_mod_k_0']] + ['total']

leloir_pathway_enz = ['EZ_' + x for x in ['UDPG4E_UDPGLUCEPIM_CPLX_mod_nad_0',
                                          'UGLT_GALACTURIDYLYLTRANS_CPLX_mod_fe2_mod_zn2_0',
                                          'LACZpp_EG12013_MONOMER',
                                          'LACZ_BETAGALACTOSID_CPLX_mod_mg2',
                                          'LCTSt3ipp_B2170_MONOMER',
                                          'LCTSt3ipp_B0070_MONOMER',
                                          'LCTSt3ipp_YDEA_MONOMER',
                                          'LCTStpp_LACY_MONOMER',
                                          'GALKr_G7096_MONOMER',
                                          'GALKr_GALACTOKIN_MONOMER_mod_mg2',
                                          ]] + ['total']
glc_pathway_enz = ['EZ_' + x for x in ['HEX1_GLUCOKIN_MONOMER_0',
                                       'GLCabcpp_ABC_18_CPLX',
                                       'GLCptspp_CPLX_157',
                                       'GLCptspp_CPLX_164',
                                       'GLCptspp_CPLX_165',
                                       'GLCt2pp_GALP_MONOMER',
                                       ]] + ['total']

lacz_pathway_enz = ['EZ_' + x for x in [
                                          'LACZpp_EG12013_MONOMER',
                                          'LACZ_BETAGALACTOSID_CPLX_mod_mg2',
                                       ]] + ['total']

leloir_pathway = ['PGMT',
                  'UDPG4E',
                  'UGLT',
                  'LACZpp',
                  'LACZ',
                  'LCTSt3ipp',
                  'LCTStpp',
                  'GALKr',
                  'GALabcpp',
                  'GALt2pp',
                  ] + ['total']
glc_pathway = ['HEX1',
               'GLCabcpp',
               'GLCptspp',
               'GLCt2pp',
               ] + ['total']


groups = {'fluxes':fluxes,
          # 'glc_enzymes':glc_enzymes,
          # 'lcts_enzymes':lcts_enzymes,
          # 'gal_enzymes':gal_enzymes,
          # 'ac_enzymes':ac_enzymes,
          'glc_fluxes':glc_fluxes,
          'lcts_fluxes':lcts_fluxes,
          'gal_fluxes':gal_fluxes,
          'ac_fluxes':ac_fluxes,
          'species':species,
          'glc_pathway_enz':glc_pathway_enz,
          'leloir_pathway_enz':leloir_pathway_enz,
          'lacz_pathway_enz':lacz_pathway_enz,
          'leloir_pathway':leloir_pathway,
          'glc_pathway':glc_pathway,
          # 'ACKr_checkpoint':ackr_checkpoint,
          # 'PYK_checkpoint':pyk_checkpoint,
          # 'PFK_checkpoint':['EZ_PFK'],
          # 'PPC_checkpoint':['EZ_PPC'],
          # 'ICD_checkpoint':['EZ_ICDHyr'],
          # 'ACS_checkpoint':['EZ_ACS_ACS_MONOMER_0'],
          }

def rescale_pct(time_data):
    time_data[time_data.index.str.startswith('EZ_')] /= time_data.loc['IV_prot_ggdw']
    time_data[time_data.index.str.startswith('MR_')] /= time_data.loc['IV_mrna_ggdw']

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

    # rescale_pct(time_data)

    summarize_model(None, time_data, groups, output_path=output_path,
                    model_tag=model_tag, backend=args.backend)
