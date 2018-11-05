#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#       E. Coli. Debug Model
#5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Values reported are roughly estimated for debugging the model. These are not
# actual reviewed values.


from os.path import join as pjoin

import cobra

import pandas as pd

import numpy as np
import logging


from pytfa.io import load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data
from pytfa.optim.relaxation import relax_dgo

from etfl.core import Enzyme, Ribosome, RNAPolymerase, ThermoMEModel, MEModel
from etfl.core.mrna import mRNA

from etfl.io.json import save_json_model

from collections import defaultdict

import re

data_dir = '../organism_data/info_ecoli'

vanilla_model = cobra.io.load_json_model('iJO1366_with_xrefs.json')

solver = 'optlang-gurobi'
# solver = 'optlang-cplex'

# Relax ATPM
# vanilla_model.reactions.ATPM.lower_bound = 0

# Load Thermo data
thermo_data = load_thermoDB('../../pytfa/data/thermo_data.thermodb')
lexicon = read_lexicon('thermo_data/iJO1366_lexicon.csv')
# lexicon = curate_lexicon(read_lexicon('thermo_data/iJO1366_lexicon.csv'))
compartment_data = read_compartment_data('thermo_data/iJO1366_compartment_data.json')

# McCloskey2014 values
glc_uptake = 7.54
glc_uptake_std = 0.56
observed_growth = 0.61 - 0.02

vanilla_model.reactions.EX_glc__D_e.lower_bound = -1*glc_uptake - glc_uptake_std
vanilla_model.reactions.EX_glc__D_e.upper_bound = -1*glc_uptake + glc_uptake_std

growth_reaction_id = 'BIOMASS_Ec_iJO1366_WT_53p95M'


#------------------------------------------------------------
# Initialisation
#------------------------------------------------------------
vanilla_model.objective = growth_reaction_id
fba_sol = vanilla_model.optimize()
mu_0 = fba_sol.f
mu_range = [0, 4]
n_mu_bins = 256

# Initialize the cobra_model
ecoli = ThermoMEModel(thermo_data,model = vanilla_model,
# ecoli = MEModel(model = vanilla_model,
                growth_reaction = growth_reaction_id,
                mu_range = mu_range,
                n_mu_bins = n_mu_bins,
                max_enzyme_concentration = 1,
                scaling = 1e3
                # scaling = 1e5
                )
ecoli.name = 'tutorial'
ecoli.logger.setLevel(logging.WARNING)

# Set the solver
ecoli.solver = solver
ecoli.solver.configuration.verbosity = 1
ecoli.solver.configuration.tolerances.feasibility = 1e-9
if solver == 'optlang_gurobi':
    ecoli.solver.problem.Params.NumericFocus = 3
ecoli.solver.configuration.presolve = True


# ------------------------------------------------------------
# Thermo
# ------------------------------------------------------------

def curate_lexicon(lexicon):
    ix = pd.Series(lexicon.index)
    ix = ix.apply(lambda s: str.replace(s,'-','__'))
    ix = ix.apply(lambda s: '_'+s if s[0].isdigit() else s)
    lexicon.index = ix
    return lexicon

lexicon = curate_lexicon(read_lexicon('thermo_data/iJO1366_lexicon.csv'))

# Annotate the cobra_model
annotate_from_lexicon(ecoli, lexicon)
apply_compartment_data(ecoli, compartment_data)

# TFA conversion
ecoli.prepare()
ecoli.convert()#add_displacement = True)


#------------------------------------------------------------
# Data
#------------------------7.54-------------------------------------

# Growth-related abundances

neidhardt_data = pd.read_excel(pjoin(data_dir,'neidhardt_tab2.xlsx'),
                               skiprows=range(0,6),
                               skip_footer=22)
mu_cols = ['mu=0.6','mu=1.0','mu=1.5','mu=2.0','mu=2.5']
neidhardt_data.columns = ['parameter','symbol','units',*mu_cols,
                          'observed_parameters','footnote']
neidhardt_data.set_index('symbol', inplace=True)

Pc = neidhardt_data.loc['Pc (μg)'][mu_cols] # μg/10^9 cells
Rc = neidhardt_data.loc['Rc (μg)'][mu_cols] # μg/10^9 cells
Mc = neidhardt_data.loc['Mc (μg)'][mu_cols] # μg dry weight/10^9 cells

neidhardt_prel = (Pc/Mc).astype(float)
neidhardt_rrel = (Rc/Mc).astype(float)
neidhardt_mu = pd.Series(Pc.index.str.replace('mu=','')).astype(float)

#------------------------------------------------------------
# Expression
#------------------------------------------------------------

# Data
# Sequences from KEGG
nt_sequences = pd.read_csv('iJO1366_nt_seq_kegg.csv',
                           index_col = 0,
                           header = None)[1]
# iJO kcat info from:
# Davidi, Dan, et al.
# "Global characterization of in vivo enzyme catalytic rates and their correspondence to in vitro kcat measurements."
# Proceedings of the National Academy of Sciences 113.12 (2016): 3401-3406.
kcat_info_milo = pd.read_excel('../models/pnas.1514240113.sd01.xlsx',
                               sheet_name='kcat 1s',
                               header=1,
                               )
kcat_info_aggregated    = pd.read_csv(pjoin(data_dir,'aggregated_kcats.csv'),
                                      index_col = 0)
ec_info_ecocyc          = pd.read_csv(pjoin(data_dir,'complex2ec.csv'),
                                      index_col = 0)
composition_info_ecocyc = pd.read_csv(pjoin(data_dir,'complex2genes.csv'),
                                      index_col = 0)
reaction2complexes_info_obrien = pd.read_excel(
    pjoin(data_dir, 'obrien2013_SI_tab10.xlsx'), index_col=0, usecols=[0, 1])
complexes2peptides_info_obrien = pd.read_excel(
    pjoin(data_dir, 'obrien2013_SI_tab1.xlsx'), index_col=0, usecols=[0, 1])

gene_names = pd.read_csv(pjoin(data_dir,'gene2bname.txt'), delimiter='\t',
                         index_col=0)


# mRNA degardation rates from
# Bernstein et al. (2002) Proc. Natl. Acad. Sci. USA, 10.1073/pnas.112318199
# "Global analysis of mRNA decay and abundance in Escherichia coli at single-gene resolution using two-color fluorescent DNA microarrays"
bernstein_ecoli_deg_rates = pd.read_excel(
    pjoin(data_dir,'bernstein_2002_mrna_deg.xls'),
    skiprows=range(8),
    index_col=0)


# Bionumbers
# ID        104876
# Property  Amino acid composition of the proteins from E. coli cell supernatant
# http://bionumbers.hms.harvard.edu/bionumber.aspx?id=104876
# Unit: per 100 moles aas
aa_ratios = {
    'K':9.01/100,  # lys__L_c
    'H':1.91/100,  # his__L_c
    'R':7.30/100,  # arg__L_c
    'D':8.30/100,  # asp__L_c
    'E':10.08/100, # glu__L_c
    'G':8.18/100,  # gly_c
    'A':10.98/100, # ala__L_c
    'V':9.63/100,  # val__L_c
    'L':7.40/100,  # leu__L_c
    'I':5.51/100,  # ile__L_c
    'T':5.22/100,  # thr__L_c
    'S':4.38/100,  # ser__L_c
    'P':3.67/100,  # pro__L_c
    'Y':1.78/100,  # tyr__L_c
    'F':3.03/100,  # phe__L_c
    'W':0.69/100,  # trp__L_c
    'M':2.40/100,  # met__L_c
    'C':0.53/100,  # cys__L_c
}

aa_dict = {
    'A':'ala__L_c',
    'R':'arg__L_c',
    'N':'asn__L_c',
    'D':'asp__L_c',
    #'B':'asx',
    'C':'cys__L_c',
    'E':'glu__L_c',
    'Q':'gln__L_c',
    #'Z':'glx',
    'G':'gly_c',
    'H':'his__L_c',
    'I':'ile__L_c',
    'L':'leu__L_c',
    'K':'lys__L_c',
    'M':'met__L_c',
    'F':'phe__L_c',
    'P':'pro__L_c',
    'S':'ser__L_c',
    'T':'thr__L_c',
    'U':'selcys__L_c', # Let's just assume selenocysteine is cysteine
    'W':'trp__L_c',
    'Y':'tyr__L_c',
    'V':'val__L_c',
    }

# TODO: Get the actual number
nt_ratios = {'u':0.25,
             'a':0.25,
             'g':0.25,
             'c':0.25,
             }

rna_nucleotides = {'u':'ura_c',
                   'g':'gua_c',
                   'a':'ade_c',
                   'c':'csn_c'}

coupling_dict = dict()

# Add cystein -> selenocystein transformation for convenience
selcys = cobra.Metabolite(id='selcys__L_c', compartment = 'c')
selcys_rxn = cobra.Reaction(id='PSEUDO_selenocystein_synthase',
                            name='PSEUDO Selenocystein_Synthase')
selcys_rxn.add_metabolites({ecoli.metabolites.cys__L_c:-1,
                            selcys:+1})
ecoli.add_reactions([selcys_rxn])

# Prot degradation
# http://www.jbc.org/content/246/22/6956.full.pdf
# The total amount of enzyme undergoing degrada- tion (2 to 7%) was the same
# during growth and during various kinds of starvation.
kdeg_low, kdeg_up = 0.02, 0.07
kdeg_enz = (kdeg_low + kdeg_up)/2

# From :
# http://book.bionumbers.org/how-fast-do-rnas-and-proteins-degrade/
# Figure 1: Measured half lives of mRNAs in E. coli, budding yeast and mouse NIH3T3 fibroblasts.
# (A, adapted from J. A. Bernstein et al., Proc. Natl Acad. Sci. USA 99:9697, 2002;
#  B, adapted from Y. Wang et al., Proc. Natl Acad. Sci. USA 99:5860, 2002;
#  C. adapted from B. Schwanhausser, Nature, 473:337, 2013).
# -------
# Mean half life of mrna is 5 minutes in ecoli
# tau = t_0.5 /ln(2)
# kdeg = 1-exp(1hr/tau)
kdeg_mrna = 1-np.exp(-60*np.log(2)/5)

# Average mrna length from Bionumber 100023
# http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=100023&ver=3
mrna_length_avg = 370

# Average peptide length
peptide_length_avg = int(np.round(mrna_length_avg/3))


# Generate a coupling dict
def is_gpr(s):
    return bool(s) and s != '[]'


# Milo kcats
#############
for x in ecoli.reactions:
    data = kcat_info_milo[kcat_info_milo['reaction (model name)'] == x.id]
    data_reverse = kcat_info_milo[kcat_info_milo['reaction (model name)'] == x.id + '_reverse']

    kcat_fwd = 0
    kcat_bwd = 0
    composition = {}

    if len(data)>0:
        candidate_complexes = data.iloc[0]
        kcat_fwd = candidate_complexes['kcat per active site [1/s]'] \
                   * candidate_complexes['catalytic sites per complex'] \
                   * 3600 # s/h
        composition = {candidate_complexes['bnumber']:candidate_complexes['polypeptides per complex']}

    if len(data_reverse)>0:
        this_data_reverse = data_reverse.iloc[0]
        kcat_bwd = this_data_reverse['kcat per active site [1/s]'] \
                * this_data_reverse['catalytic sites per complex'] \
                * 3600 # s/h
        composition = {
            this_data_reverse['bnumber']:this_data_reverse['polypeptides per complex']} \
                if not composition else composition

    if kcat_fwd == 0 and kcat_bwd == 0:
        continue

    if kcat_bwd == 0:
        kcat_bwd = kcat_fwd

    if kcat_fwd == 0:
        kcat_fwd = kcat_bwd

    #FIXME several polypeptides per complex ??

    new_enzyme = Enzyme(x.id,
                        kcat_fwd=kcat_fwd,
                        kcat_bwd=kcat_bwd,
                        kdeg=kdeg_enz)

    new_enzyme.composition = composition

    coupling_dict[x.id] = [new_enzyme]


aggregated_coupling_dict = defaultdict(list)


# Aggregated kcats
##################

def ec2ecocyc(ec_number):
    if not isinstance(ec_number, list):
        return ec_info_ecocyc[ec_info_ecocyc['ec'] == ec_number]
    else:
        return ec_info_ecocyc[ec_info_ecocyc['ec'].isin(ec_number)]

def score_against_genes(putative_genes, reaction_genes):
    score = 0

    putative_genes_list = putative_genes.split('" // "')
    reaction_gene_list  = [x.name for x in reaction_genes]
    s1 = len(set(putative_genes_list).intersection(reaction_gene_list))
    s2 = len(set(putative_genes_list).difference  (reaction_gene_list))
    s3 = len(set(reaction_gene_list) .difference  (putative_genes_list))

    score = 2*s1-s2-s3
    # print(putative_genes_list, reaction_gene_list, score)
    return score

def match_ec_genes_ecocyc(ecocyc, genes, threshold=0.5):
    this_data = composition_info_ecocyc[composition_info_ecocyc['complex'].isin(ecocyc)]
    scores = this_data['putative_genes'].apply(score_against_genes, args=[genes])
    selectable = this_data[scores>len(genes)*threshold]
    if len(selectable) == 0:
        return None, scores
    else:
        return selectable[scores == scores.max()].iloc[0], scores

def ecocyc2composition(ecocyc):
    ecocyc_comp = composition_info_ecocyc[composition_info_ecocyc['complex'] == ecocyc]
    # We do a left join to get the bnumbers that are used in the model
    composition_data = ecocyc_comp.merge(gene_names, right_index=True,
                                  how='left', left_on='gene')
    composition = defaultdict(int)
    for e,row in ecocyc_comp.iterrows():
        this_gene = row['gene']
        try:
            this_b_number = gene_names.loc[this_gene]['b#']
            composition[this_b_number] += composition_data['coeffs'].iloc[0]
        except:
            ecoli.logger.warning('Could not find gene associated to {}'
                                 .format(row['obj_ids']))
            ecoli.logger.info(ecocyc_comp)

    return composition

comp_regex = re.compile(r'(b[0-9]{4})\((\d?)\)')

def complex2composition(complex_name):
    # Silence modifications
    if '_mod_' in complex_name:
        complex_name = complex_name[0:complex_name.index('_mod_')]

    composition_string = complexes2peptides_info_obrien.loc[complex_name,'Gene composition']
    composition_dict = {}
    groups = comp_regex.findall(composition_string)
    for peptide, stoich in groups:
        if stoich == '':
            stoich = 1
        else:
            stoich = int(stoich)
        composition_dict[peptide] = stoich

    return composition_dict

def ec2kcat(ec_number):
    try:
        return kcat_info_aggregated['kcat'].loc[ec_number].max() * 3600  # s/h
    except KeyError:
        return None


for x in ecoli.reactions:
    if x.id in coupling_dict:
        # We already have info
        continue

    if 'ec_numbers' not in x.notes or x.notes['ec_numbers'] == ['nan']:
        # There is nothing we can do
        continue

    reaction_ecs = x.notes['ec_numbers']

    try:
        complex_names = reaction2complexes_info_obrien.loc[x.id,'Enzymes'].split(' OR ')
    except KeyError:
        continue

    for e,this_complex_name in enumerate(complex_names):

        # Start with this:
        composition = complex2composition(this_complex_name)
        if not composition:
            # Skip this one
            continue

        this_ec = x.notes['ec_numbers'][0]
        kcat = ec2kcat(this_ec)

        if kcat is None:
            continue

        new_enzyme = Enzyme('{}_{}'.format(x.id,this_complex_name),
                            name='{}_{}: {}'.format(x.id, e, this_complex_name),
                            kcat=kcat,
                            kdeg=kdeg_enz)

        new_enzyme.composition = composition

        new_enzyme.notes['EC'] = this_ec

        aggregated_coupling_dict[x.id].append(new_enzyme)


# 1/0
coupling_dict.update(aggregated_coupling_dict)
mrna_dict = dict()

# Generate a mRNA dict

for x in nt_sequences.index:
    try:
        the_gene = ecoli.genes.get_by_id(x)
    except KeyError:
        ecoli.add_genes([cobra.Gene(id=x)])
        the_gene = ecoli.genes.get_by_id(x)

    # Try to get half life from Bernstein et al.
    try:
        t_half = bernstein_ecoli_deg_rates.loc[x.upper()]['medium, min.1'] #M9 medium
        # Mean half life of mrna is 5 minutes in ecoli
        # tau = t_0.5 /ln(2)
        # kdeg = 1-exp(1hr/tau)
        this_kdeg_mrna = 1 - np.exp(-60 * np.log(2) / t_half)
    except KeyError:
        this_kdeg_mrna = kdeg_mrna # Average value of 5 mins

    if np.isnan(this_kdeg_mrna):
        this_kdeg_mrna = kdeg_mrna # Average value of 5 mins

    new_mrna = mRNA(x,`
                    kdeg = this_kdeg_mrna,
                    gene_id = the_gene.id)
    mrna_dict[x] = new_mrna

#[ ribosomes and RNAP ]#
rpeptide_genes = pd.read_csv(pjoin(data_dir,'ribosomal_proteins_ecoli.tsv'),
                             delimiter='\t',
                             header=None)[0]
rpeptide_genes = rpeptide_genes.str.split(':').apply(lambda x:x[1])

rib = Ribosome(id='rib', name='Ribosome', kribo=12 * 3600, kdeg=0.001)
rnap = RNAPolymerase(id='rnap',
                     name='RNA Polymerase',
                     ktrans = 1000*3600,
                     kdeg = 0.2)


##########################
##    MODEL CREATION    ##
##########################

ecoli.add_rnap(rnap)
ecoli.add_ribosome(rib, free_ratio=0.2)
ecoli.add_nucleotide_sequences(nt_sequences)
ecoli.add_mrnas(mrna_dict.values())
# RNAP
# b3295: K03040 DNA-directed RNA polymerase subunit alpha [EC:2.7.7.6] | (RefSeq) rpoA; RNA polymerase, alpha subunit
# b3649: K03060 DNA-directed RNA polymerase subunit omega [EC:2.7.7.6] | (RefSeq) rpoZ; RNA polymerase, omega subunit
# b3987: K03043 DNA-directed RNA polymerase subunit beta [EC:2.7.7.6] | (RefSeq) rpoB; RNA polymerase, beta subunit
# b3988: K03046 DNA-directed RNA polymerase subunit beta' [EC:2.7.7.6] | (RefSeq) rpoC; RNA polymerase, beta prime subunit
# Ribosome
# rRNA:
# b3851: K01977 16S ribosomal RNA | (RefSeq) rrsA; 16S ribosomal RNA of rrnA operon
# b3854: K01980 23S ribosomal RNA | (RefSeq) rrlA; 23S ribosomal RNA of rrnA operon
# b3855: K01985 5S ribosomal RNA | (RefSeq) rrfA; 5S ribosomal RNA of rrnA operon
# rPeptides:
# See file ribosomal_proteins_ecoli.tsv
ecoli.build_expression(aa_dict = aa_dict,
                       rna_nucleotides= rna_nucleotides,
                       atp='atp_c',
                       amp='amp_c',
                       gtp='gtp_c',
                       gdp='gdp_c',
                       ppi='ppi_c',
                       h2o='h2o_c',
                       h='h_c',
                       rnap_genes = ['b3295','b3649','b3987','b3988'],
                       rrna_genes = ['b3851','b3854','b3855'],
                       rprot_genes= rpeptide_genes
                       )
ecoli.add_enzymatic_coupling(coupling_dict)
ecoli.populate_expression()
ecoli.add_degradation()
ecoli.add_interpolation_variables()
ecoli.add_dummies(nt_ratios=nt_ratios,
                  mrna_kdeg=kdeg_mrna,
                  mrna_length=mrna_length_avg,
                  aa_ratios=aa_ratios,
                  enzyme_kdeg=kdeg_enz,
                  peptide_length=peptide_length_avg)
ecoli.add_protein_mass_requirement(neidhardt_mu, neidhardt_prel)
ecoli.add_rna_mass_requirement(neidhardt_mu, neidhardt_rrel)


ecoli.print_info()
ecoli.optimize()
# ecoli.growth_reaction.lower_bound = observed_growth
# relaxed_model, slack_model, relax_table = relax_dgo(ecoli)
#
# solution = relaxed_model.optimize()
# print('Growth               : {}'.format(relaxed_model.solution.f))
# print(' - Ribosomes produced: {}'.format(relaxed_model.solution.x_dict.EZ_rib))
# print(' - RNAP produced: {}'.format(relaxed_model.solution.x_dict.EZ_rnap))
#
# filepath = 'models/iJO1366_t_{}_{}_bins_20180503.json'.format(
#     len(relaxed_model.enzymes), relaxed_model.n_mu_bins)
# save_json_model(relaxed_model, filepath)
#
#
