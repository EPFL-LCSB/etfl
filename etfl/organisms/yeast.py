from os.path import join as pjoin
import os

import cobra

import pandas as pd

import numpy as np
from numbers import Number

from pytfa.io import load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data
from pytfa.io.base import import_matlab_model

from .ecoli_utils import infer_enzyme_from_gpr

from ..core.expression import is_me_compatible
from ..core import Enzyme, Ribosome, RNAPolymerase, ThermoMEModel, MEModel
from ..core.rna import mRNA

from ..optim.constraints import ConstantAllocation
from ..optim.variables import EnzymeVariable
from pytfa.optim.utils import symbol_sum

from collections import defaultdict
from warnings import warn
import re

from .yeast_utils import compositions_from_string, special_coupling, find_bbb_ratio, is_transport


partial_saturation = False
if partial_saturation:
    sigma = 0.51
else:
    sigma = 1


def clean_string(s):

    s = s.replace('-','_')

    # Remove invalid characters
    s = re.sub('[^0-9a-zA-Z_]', '', s)

    # Remove leading characters until we find a letter or underscore
    s = re.sub('^[^a-zA-Z_]+', '', s)

    return s


#data_dir = '/src/yetfl/data'
data_dir = '../organism_data/info_yeast'


#########################
###     BASE MODEL    ###
#########################

def get_model(solver):
    vanilla_model = cobra.io.load_matlab_model(pjoin(data_dir,'Y8_3_4_mod_curated.mat'))
#    vanilla_model = cobra.io.load_matlab_model('../project_debug/allocation/allocation_mod.mat')
    vanilla_model.slim_optimize()
    vanilla_model.solver = solver
    vanilla_model.slim_optimize()
    


    def sanitize_varnames(model):
        for met in model.metabolites:
            if met.id[0].isdigit():
                met.id = '_' + met.id
        for rxn in model.reactions:
            if rxn.id[0].isdigit():
                rxn.id = '_' + rxn.id
        model.repair()

        return model

    # Add cystein -> selenocystein transformation for convenience
#    selcys = cobra.Metabolite(id='selcys_L_c', compartment='c', formula='C3H7NO2Se')
#    selcys_rxn = cobra.Reaction(id='PSEUDO_selenocystein_synthase',
#                                name='PSEUDO Selenocystein_Synthase')
#    selcys_rxn.add_metabolites(
#        {vanilla_model.metabolites.get_by_id('s_0981_c'): -1, selcys: +1})
#    vanilla_model.add_reactions([selcys_rxn])

    sanitize_varnames(vanilla_model)
    vanilla_model.slim_optimize()
    return vanilla_model

def get_comp_from_id(met_id):
    splt_1 = met_id.split('[')
    str_1 = splt_1[1]
    splt_2 = str_1.split(']')
    str_2 = splt_2[0]
    # a dictionary to convert compartment id from Y8 to Y7
    compID={'ce' : 'q',
            'er' : 'r',
            'p' : 'x' ,
            'vm' : 'w',
            'mm' : 'o',
            'gm' : 'j',
            'lp' : 'l',
            'erm': 't',}

    return compID[str_2]

# ------------------------------------------------------------
# Thermo
# ------------------------------------------------------------

def get_thermo_data():
    # Load Thermo data
    thermo_data = load_thermoDB(pjoin(data_dir,'thermo/thermo_data.thermodb'))
    lexicon = read_lexicon(pjoin(data_dir,'thermo/yeast_lexicon.csv'))
    compartment_data = read_compartment_data(pjoin(data_dir,'thermo/yeast_compartment_data.json'))


    def curate_lexicon(lexicon):
        ix = pd.Series(lexicon.index)
        ix = ix.apply(lambda s: str.replace(s,'-','__'))
        ix = ix.apply(lambda s: '_'+s if s[0].isdigit() else s)
        lexicon.index = ix
        return lexicon

    lexicon = curate_lexicon(lexicon)

    return thermo_data, lexicon, compartment_data

#------------------------------------------------------------
# Data
#------------------------7.54--------------------------------

# Essentials
def get_essentials():
    return dict(atp='s_0434_c',
                adp='s_0394_c',
                amp='s_0423_c',
                gtp='s_0785_c',
                gdp='s_0739_c',
                pi='s_1322_c',
                ppi='s_0633_c',
                h2o='s_0803_c',
                h='s_0793_c')
# growth-related abundances    
def get_mass_ratios():
    # H. C. Lange, J. J. Heijnen, 2001, Statistical Reconciliation of the
    # Elemental and Molecular Biomass Composition of Saccharomyces cerevisiae
    ratios = pd.read_csv(pjoin(data_dir,'mass_ratios.csv'),header=None, index_col=0)
    mu_cols = ['mu=0.022','mu=0.052','mu=0.087','mu=0.107','mu=0.126', \
               'mu=0.158','mu=0.211','mu=0.3','mu=0.35', \
                   'mu=0.4']
    ratios.columns = mu_cols
    
    mu = ratios.loc['mu']
    rna = ratios.loc['RNA']
    prot = ratios.loc['Protein']
    dna = ratios.loc['DNA']
    lipid = ratios.loc['Lipid']
    carbohydrate = ratios.loc['Carbohydrate']
    ion = ratios.loc['Ion']
    water = ratios.loc['Water']
    pi = ratios.loc['Pi']
    total = ratios.loc['Sum']
    
    new_total = (total - pi- water).astype(float) # remove water and pi from the composition
    
    prot = (prot/new_total).astype(float)
    rna = (rna/new_total).astype(float)
    dna = (dna/new_total).astype(float)   
    carbohydrate = (carbohydrate/new_total).astype(float)   
    lipid = (lipid/new_total).astype(float)
    ion = (ion/new_total).astype(float)   

    return mu, rna, prot, dna, lipid, carbohydrate, ion


#------------------------------------------------------------
# Expression
#------------------------------------------------------------

# Data
# Sequences from KEGG
nt_sequences = pd.read_csv(pjoin(data_dir,'Y8_nt_seq_kegg.csv'),
                           index_col = 0,
                           header = None)[1]

def get_nt_sequences():
    return nt_sequences

# Sequences from KEGG
aa_sequences = pd.read_csv(pjoin(data_dir,'Y8_aa_seq_kegg.csv'),
                           index_col = 0,
                           header = None)[1]

def get_aa_sequences():
    return aa_sequences

# Bionumber
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=102126&ver=2&trm=Saccharomyces+cerevisiae+gc+ratio&org=
# ID	102126
# Property	GC content 
# Organism	Budding yeast Saccharomyces cerevisiae
# Value	38.3
# Units	%
# Reference	Wood et al. The genome sequence of Schizosaccharomyces pombe. 2002. Nature 415, 871-880. P.872 right column 4th paragraph
# Assumed GC ratio of 0.383
gc_ratio = 0.383 # %
#http://www.ensembl.org/Saccharomyces_cerevisiae/Location/Genome
chromosome_len = 12157105 # bp

def get_yeast_gen_stats():
    return chromosome_len, gc_ratio

def get_ratios():
    # Bionumbers
    # ID        112800
    # Property  Amino acid composition of the proteins as measured
    # https://bionumbers.hms.harvard.edu/bionumber.aspx?id=112800&ver=0&trm=Saccharomyces+cerevisiae+amino+acid+ratio&org=
    # Unit: mole %
    aa_ratios = {
        'K':6.57/100,  # lys__L_c
        'H':1.93/100,  # his__L_c
        'R':3.86/100,  # arg__L_c
        'D':15.48/100, # asp__L_c
        'E':10.08/100, # glu__L_c
        'G':8.89/100,  # gly_c
        'A':9.77/100,  # ala__L_c
        'V':7.33/100,  # val__L_c
        'L':8.01/100,  # leu__L_c
        'I':5.89/100,  # ile__L_c
        'T':5.57/100,  # thr__L_c
        'S':5.33/100,  # ser__L_c
        'P':4.22/100,  # pro__L_c
        'Y':1.96/100,  # tyr__L_c
        'F':3.76/100,  # phe__L_c
        'W':0.65/100,  # trp__L_c
        'M':1.14/100,  # met__L_c
        'C':0.14/100,  # cys__L_c
        #'O':0.24/100,  # orn__L_c
    }


    nt_ratios = {'u':0.5*(1-gc_ratio),
                 'a':0.5*(1-gc_ratio),
                 'g':0.5*gc_ratio,
                 'c':0.5*gc_ratio,
                 }

    return nt_ratios,aa_ratios


def get_monomers_dict():

    aa_dict = {'A': 's_0955_c',
               'R': 's_0965_c',
               'N': 's_0969_c',
               'D': 's_0973_c',
               # 'B':'asx',
               'C': 's_0981_c',
               'E': 's_0991_c',
               'Q': 's_0999_c',
               # 'Z':'glx',
               'G': 's_1003_c',
               'H': 's_1006_c',
               'I': 's_1016_c',
               'L': 's_1021_c',
               'K': 's_1025_c',
               'M': 's_1029_c',
               'F': 's_1032_c',
               'P': 's_1035_c',
               'S': 's_1039_c',
               'T': 's_1045_c',
#               'U': 'selcys_L_c',
               'W': 's_1048_c',
               'Y': 's_1051_c',
               'V': 's_1056_c', }

    rna_nucleotides = {
                'u': 's_1559_c',
                'g': 's_0785_c',
                'a': 's_0434_c',
                'c': 's_0539_c'}

    rna_nucleotides_mp = {
                'u': 's_1545_c',
                'g': 's_0782_c',
                'a': 's_0423_c',
                'c': 's_0526_c'}

    dna_nucleotides = {
                't': 's_0650_c',
                'g': 's_0617_c',
                'a': 's_0586_c',
                'c': 's_0590_c'}


    return aa_dict, rna_nucleotides, rna_nucleotides_mp, dna_nucleotides


def remove_from_biomass_equation(model, nt_dict, aa_dict, atp_id, adp_id,
                                 h2o_id, h_id, pi_id):

    # According to discussions, should only remove GAM

    mets_to_rm = dict()

    old_total_stoich = abs(sum([x for x in
                                model.growth_reaction.metabolites.values() if
                                x<0]))

    expression_mets = list(nt_dict.values()) + list(aa_dict.values())

    for m,stoich in model.growth_reaction.metabolites.items():
        if m.id in expression_mets:
            mets_to_rm[m] = -1*stoich

    model.growth_reaction.add_metabolites(mets_to_rm)

    # We need to add back the ATP from the GAM that we just removed but should
    # still be taken into account. Indeed, nADP =/= nATP in the original
    # biomass reaction:
    # -54.119975 atp_c + .... --> 53.95 adp_c

    atp = model.metabolites.get_by_id(atp_id)
    adp = model.metabolites.get_by_id(adp_id)
    pi  = model.metabolites.get_by_id(pi_id)
    h2o = model.metabolites.get_by_id(h2o_id)
    h = model.metabolites.get_by_id(h_id)
    atp_recovery = model.growth_reaction.metabolites[adp]
    model.growth_reaction.add_metabolites({atp:-1*atp_recovery})


# Prot degradation
# Bionumbers
# ID        115093
# Property  Average and median half life of protein
# Organism  Saccharomyces cerevisiae
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=115093&ver=1&trm=protein+half+life+yeast&org=
# Unit: hour
# Value: 8.8
# tau = t_0.5 /ln(2)
# kdeg = ln(2)/t_half * (conversion to hours)
kdeg_enz = np.log(2)/8.8

# mRNA degradation
# Bionumbers
# ID        105511
# Property  Mean mRNA half life 
# Organism  Saccharomyces cerevisiae
# https://https://bionumbers.hms.harvard.edu/bionumber.aspx?id=105511&ver=7&trm=Saccharomyces+cerevisiae+mean+half+life+of+mrna&org=
# Unit: min
# Value: 23
# tau = t_0.5 /ln(2)
# kdeg = ln(2)/t_half * (conversion to hours)
kdeg_mrna = np.log(2)/(23/60)

# Average mrna length from Bionumber 107678
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=107678&ver=3&trm=Saccharomyces+cerevisiae+length+of+chromosomes&org=
# mrna_length_avg = 1250
# Unit: nucleotides
mrna_length_avg = 1250

# Average peptide length
peptide_length_avg = int(np.round(mrna_length_avg/3))


def get_mrna_metrics():
    return kdeg_mrna, mrna_length_avg

def get_enz_metrics():
    return kdeg_enz, peptide_length_avg

# Generate a coupling dict
def is_gpr(s):
    return bool(s) and s != '[]'



# Y8 kcat info from:
# GECKO
# Improving the phenotype predictions of a yeast genome-scale metabolic model
# by incorporating enzymatic constraints
# kmax_info_gecko = pd.read_excel(pjoin(data_dir,'pnas.1514240113.sd01.xlsx'),
#                                sheet_name='kmax 1s',
#                                header=2,
#                                )
#kcat_info_aggregated    = pd.read_excel(pjoin(data_dir,'aggregated_kcats.xlsx'),
#                                      sheet_name='kmax')

# data for each related protein abundance from PaxDB obtained from gecko
prot_abundance = pd.read_csv(pjoin(data_dir,'abundance_table.csv'),header=0)

# data for each protein MWs (including enzymes) from swissprot obtained from gecko
mws = pd.read_csv(pjoin(data_dir,'all_peptides_mw.csv'),header=0)

ec_info_yeastcyc          = pd.read_excel(pjoin(data_dir,'complex2kcat.xlsx'),
                                      sheet_name = 'complex2kcat')
# ec_info_yeastcyc          = pd.read_excel(pjoin(data_dir,'enzyme2kcat.xlsx'),
#                                       sheet_name = 'enzyme2kcat')
# Gecko: the observed saturation rate for the enzymes is 0.51
ec_info_yeastcyc['kcat'] *= sigma
predicted_enzyme_ids = ec_info_yeastcyc['complex'] # this is only important for enzymes with multiple kcats
# Gecko kcat manual curation
# kcat_manual = pd.read_csv(pjoin(data_dir,'manual_curation.csv'),header=0)
# kcat_manual['kcat'] *= sigma

gene_list =  pd.read_excel(pjoin(data_dir,'trans_list.xlsx'),
                                      sheet_name = 'Genes')
pep_list =  pd.read_excel(pjoin(data_dir,'trans_list.xlsx'),
                                      sheet_name = 'Peptides')

Kcat_big = 1e+9 # a big number similar to infinity
Kcat_average = sigma*70.9*3600 # an average kcat for enzymes with missing kcat

# I wanted to integrate backward kcats from gecko, but at the end only for 2 reactions they were different
bkwd_kcats = {'YLR354C_enzyme' : 75*3600, # transaldolase
              'YLR304C_enzyme' : 143.3026*3600 # cis-aconitate(3-) to isocitrate 
              }
for k in bkwd_kcats.keys():
    bkwd_kcats[k] *= sigma

# yeastcyc kcats
################

TUPLE_TOKEN = 'TUPLE // '
SEP_TOKEN = ' // '
BROKEN_SEP_TOKEN = '"{}"'.format(SEP_TOKEN) # Looks like '" // "'

def clean_components(component_string):
    """
    Transforms an entry like :
    "TUPLE // CPLX-7491 // 2 // TUPLE // HS07692-MONOMER // 6"
    to
    {'CPLX-7491':'2',
     'HS07692-MONOMER':'6'}

    :param component_string:
    :type component_string: str
    :return components_dict:
    :type components_dict: dict(str:str)
    """
    # The first one will always be empty
    components = component_string.split(TUPLE_TOKEN)[1:]
    components_dict = dict([x.split(SEP_TOKEN)[:2] for x in components])
    return components_dict

def align_genes_products(genes, products):
    """
    the genes are in a string that looks like
    'SUCLG1" // "SUCLA2'
    (The double quotes are not a typo)
    the products are in a string that looks like
    'HS08877-MONOMER // ENSG00000136143-MONOMER'

    returns a dict of the shape:
    {'HS08877-MONOMER': 'SUCLG1'
     'ENSG00000136143-MONOMER': 'SUCLA2'}

    :param genes:
    :type genes: str
    :param products:
    :type products: str
    :return product_genes_dict:
    :type product_genes_dict: dict(str:str)
    """
    product_genes_dict = dict(zip(products.split(SEP_TOKEN),
                                  genes.split(SEP_TOKEN)))
    return product_genes_dict

def concatenate(df1,df2):
    '''attempts to concatenate two databases df1, df2,
    while avoiding including shared information between
    databases more than once.'''
    
    genes_dict1 = dict()
    genes_dict2 = dict()
    for _,row in df1.iterrows():
        this_id = row['prot_id']
        gene_name = row['gene_id']
        if isinstance(gene_name,np.float) and np.isnan(gene_name):
            continue
        genes_dict1[this_id]= tuple(gene_name.split(SEP_TOKEN))
        
    for _,row in df2.iterrows():
        this_id = row['prot_id']
        gene_name = row['gene_id']
        if isinstance(gene_name,np.float) and np.isnan(gene_name):
            continue
        this_gene_set = tuple(gene_name.split(SEP_TOKEN))
        score_df = pd.Series({k:jaccard_score(this_gene_set,v) for k,v in genes_dict1.items()})
        perfect_scores = score_df[abs(score_df-1) < 1e-4]
        if len(perfect_scores)==0:
            genes_dict2[this_id]= tuple(gene_name.split(SEP_TOKEN))
    DF = df2[df2['prot_id'].isin(genes_dict2.keys())]
    df = pd.concat([df1,DF], axis=0)
    return df

def _infer_database_information(df):
    enz_dict = dict()
    genes_dict = dict()

    for _,row in df.iterrows():
        this_id = row['prot_id']
        cpids = row['components_ids']
        gene_name = row['gene_id']
        if isinstance(cpids,np.float) and np.isnan(cpids):
            continue
        elif isinstance(gene_name,np.float) and np.isnan(gene_name):
            continue
        components = clean_components(cpids)
        product_genes_dict = align_genes_products(gene_name,
                                                  row['products_ids'])
        composition = {product_genes_dict[k]:int(v)
                             for k,v in components.items()
                             if k in product_genes_dict} # Complexes are not a
        if len(composition) == 0:
            warn('No composition found for enzyme {} with genes {}'
                 .format(this_id, gene_name))
            continue
                                                    # direct product of genes
        the_enz = Enzyme(   id = this_id,
                            kcat=get_average_kcat(),
                            kdeg=kdeg_enz,
                            composition=composition)
        enz_dict[this_id] = the_enz
        genes_dict[this_id]= tuple(gene_name.split(SEP_TOKEN))
        
    return enz_dict, genes_dict


def get_yeast_coupling_dict(model):
    # df1 =  pd.read_csv(pjoin(data_dir,'yeastcyc_prot_complexes_common.csv'),header=0)
    df1 = pd.read_csv(pjoin(data_dir,'Yeast_complex_portal.csv'),header=0)
    df2 =  pd.read_csv(pjoin(data_dir,'yeastcyc_prot_complexes_ids.csv'),header=0)
    df3 =  pd.read_csv(pjoin(data_dir,'simplified_composition.csv'),header=0)

    # df1.columns = ['prot_name','gene_name','components','products']
    df1.columns = ['prot_id','gene_id','components_ids','products_ids']
    df2.columns = ['prot_id','gene_id','components_ids','products_ids']
    df3.columns = ['prot_id','gene_id','components_ids','products_ids']
        
    # df = pd.concat([df1,df2], axis = 1)
    # Concatenating two databases while avoiding redundancy:
    df = concatenate(df1,df2)
    # the 3rd database is manually created, so the elements are not present in the other two
    df = pd.concat([df,df3], axis=0)
    
    enz_dict, genes_dict = _infer_database_information(df)        

    print('Size of the enzyme dictionary:', len(enz_dict))

    # match with Jaccard test reactions to the enz_dict
    final_dict = dict()
    
    
    new_name='{}_enzyme'
    tol = 1e-4
    half_score = 1/2
    third_score = 1/3
    

    # Now, we can try to expand the enzyme coverage for the reactions
    # that we have information about only one of the enzymes catalyzing
    # the reaction. For the other enzyme, we use stoichiometric data from
    # the known enzyme to infer its composition.
    def _two_enzymes_for_a_reaction(matched_enz, gene_ids):
        semi_perfect_scores = score_df[abs(score_df-half_score) < tol]
        semi_matched_enz=[enz_dict[x] for x in semi_perfect_scores.index]
        if len(semi_matched_enz) > 0:
            if len(semi_matched_enz) == 2:
                matched_enz=matched_enz+semi_matched_enz
            elif len(semi_matched_enz) == 1:
                # suppose that we have only one enzyme with known composition for this reaction
                semi_matched_enz=semi_matched_enz[0]
                matched_enz = matched_enz+[semi_matched_enz]
                gene_ids=list(set(gene_ids)-set(list(semi_matched_enz.composition.keys())))
                # Let's not add redundant enzymes
                gene_ids=list(set(gene_ids)-set([k.split('_')[0] for k,_ in enz_dict.items()]))
            # define a new enzyme for the gene with unknown composition
                for i in gene_ids:    
                    the_enz = Enzyme( id=new_name.format(i),
                                    kcat=get_average_kcat(),
                                    kdeg=kdeg_enz,
                                    composition={i:int(v) for _,v in semi_matched_enz.composition.items()})
                    enz_dict[new_name.format(i)] = the_enz
                    genes_dict[new_name.format(i)]=(i,)
                    matched_enz=matched_enz+[enz_dict[new_name.format(i)]]
        return enz_dict, genes_dict, matched_enz, gene_ids
    
    # It's possible that we are missing a gene in the GEM model,
    # while we have information about it in the yeast_cyc data:
    def _missing_gene_in_gem(matched_enz):
        if abs(top_scores-third_score) < tol :
            good_scores = score_df[abs(score_df-top_scores) < tol]
            sus_enz = [enz_dict[x] for x in good_scores.index]
            if len(sus_enz)==3:
                matched_enz=matched_enz+sus_enz
            elif len(sus_enz)==2: # if we have two enzymes for a reaction and we have information about only one of them, it's better to add none of them!
                check_score=jaccard_score(list(sus_enz[0].composition.keys()),list(sus_enz[1].composition.keys()))
                if abs(check_score-third_score) < tol : # If the unspecified gene in GEM is different between two enzymes, we prefer not to add them!
                    matched_enz=matched_enz+sus_enz
        return matched_enz
    
    def _find_monomeric_enz(this_rxn, the_gene_set):
        if len(gene_set) == 1:
            return True, gene_set
        elif len(gene_set) > 1: # excludes reactions with no assigned genes 
            # if gpr does not contain 'and' it can be considered as monomeric isozymes
            gpr = this_rxn.gene_reaction_rule
            gpr = gpr.replace('(','').replace(')','')
            if len(gpr.split(' and ')) == 1:
                return True, gene_set
            else: # possible: g1 or g2 or (g3 and g4)
                gpr = gpr.split(' or ')
                mono_genes = [x.replace(' ','') for x in gpr
                              if x.replace(' ','') in gene_set]
                if len(mono_genes) > 0:
                    return True, mono_genes
                else:
                    return False, []
        elif len(gene_set) == 0:
             return False, []       
            
    def _remove_missannotated_monomers(matched_enz):
        # if a monomeric enzyme is completely part of a complex it should be removed
        mon_id = '_enzyme' # this is the sign of monomeric enzymes
        monomers = [x for x in matched_enz if mon_id in x.id]
        complexes = [x for x in matched_enz if mon_id not in x.id]
        if len(monomers) == 0:
            return matched_enz
        else:
            for enz in monomers:
                for comp in complexes:
                    if enz.id.replace('_enzyme','') in comp.composition.keys():
                        matched_enz.remove(enz)
            return matched_enz
      
    special = ['r_0906', 'r_0916'] # list of reactions for special coupling that do not enter the function by themselves    
    for rxn in model.reactions:
        gene_ids = [x.id for x in rxn.genes]
        gene_set = gene_ids
        score_df = pd.Series({k:jaccard_score(gene_ids,v) for k,v in genes_dict.items()})
        top_scores = score_df.max()
        perfect_scores = score_df[abs(score_df-1) < tol]
           
        matched_enz = [enz_dict[x] for x in perfect_scores.index]
        if (len(perfect_scores) == 0 and top_scores > 0.49) or rxn.id in special:
            matched_enz = special_coupling(rxn, enz_dict, genes_dict)
        
        enz_dict, genes_dict, matched_enz, gene_ids = _two_enzymes_for_a_reaction(matched_enz, gene_ids)
        
        matched_enz = _missing_gene_in_gem(matched_enz)
        
        if len(matched_enz) > 0:
            matched_enz = _remove_missannotated_monomers(matched_enz) # it's possible that monomeric enzymes that partipates in complexes are missannotated here
            matched_enz = list(set(matched_enz))
            final_dict[rxn.id] = matched_enz
        # In this case, we can assume the reaction is catlayzed by a monomer,
        # but we add only monomers with a known kcat.
        else:
            det, monomeric_genes = _find_monomeric_enz(rxn, gene_set) 
            if det == True:
                enz_set = []
                for the_gene in monomeric_genes:
                    # Let's not add redundant enzymes.
                    new_enz = False
                    try:
                        enz_dict[new_name.format(the_gene)]
                    except (KeyError):
                        new_enz = True
                    if new_enz:
                        the_enz = Enzyme( id=new_name.format(the_gene),
                                            kcat=get_average_kcat(),
                                            kdeg=kdeg_enz,
                                            composition={the_gene:1})
                        enz_dict[new_name.format(the_gene)] = the_enz
                        genes_dict[new_name.format(the_gene)]=(the_gene,)
                    else:
                        the_enz = enz_dict[new_name.format(the_gene)]
                    enz_set += [the_enz]
                final_dict[rxn.id] = enz_set
        # add enzymes with multiple kcats       
        try:
            non_spec = []
            for z in final_dict[rxn.id]:
                chimeric_id = z.id + '_' + rxn.id
                multiple = predicted_enzyme_ids[predicted_enzyme_ids == chimeric_id]
                if multiple.any():
                    the_enz = Enzyme( id=chimeric_id,
                                        kcat=get_average_kcat(),
                                        kdeg=kdeg_enz,
                                        composition=z.composition)
                    non_spec += [z] 
                    final_dict[rxn.id] += [the_enz]
                    enz_dict[chimeric_id] = the_enz
            for z in non_spec:
                final_dict[rxn.id].remove(z) # need to remove non-specefic enzyme
        except KeyError:
            continue    
        # setting kcat value (until now it's the median)
        try:
            if is_transport(rxn): # for the transport rxns a big kcat is assigned just to keep them in the model but without a real catalytic constraint
                z.kcat_fwd = Kcat_big
                z.kcat_bwd = Kcat_big
            else:
                for z in final_dict[rxn.id]:
                    z.kcat_fwd=get_kcat(z.id)
                    if z.id in bkwd_kcats.keys():
                        z.kcat_bwd = bkwd_kcats[z.id]
                    else:
                        z.kcat_bwd=z.kcat_fwd
        except KeyError:
            continue
        
    print('Size of the new enzyme dictionary:', len(enz_dict))
    print('Size of the coupling dictionary:', len(final_dict))

    return  final_dict

def jaccard_score(list1, list2):
    """
    Jaccard similarity score
    The Jaccard coefficient measures similarity between finite sample sets, and
    is defined as the size of the intersection divided by the size of the union
    of the sample sets

    https://en.wikipedia.org/wiki/Jaccard_index

    :param input_gene_list:
    :param background_gene_list:
    :return:
    """
    N = len(set(list1).intersection(list2))
    D = len(set(list1).union(list2))
    if D == 0:
        return 0
    else:
        return 1.0*N/D


def get_average_kcat():
   return Kcat_average # For testing purposes


# Aggregated kcats
##################
def get_kcat(enz_id):
    kcat=complex2kcat(enz_id)
    return kcat

def complex2kcat(enz_id):

    kcat = ec_info_yeastcyc[ec_info_yeastcyc['complex'] == enz_id]['kcat']
    if not kcat.any():
        kcat = Kcat_average
    return float(kcat) 
    

def ec2complex(ec_number):
    if not isinstance(ec_number, list):
        return ec_info_yeastcyc[ec_info_yeastcyc['ec'] == ec_number]
    else:
        return ec_info_yeastcyc[ec_info_yeastcyc['ec'].isin(ec_number)]

def complex2ec(complex_name):
    if not isinstance(complex_name, list):
        return ec_info_yeastcyc[ec_info_yeastcyc['complex'] == complex_name]
    else:
        return ec_info_yeastcyc[ec_info_yeastcyc['complex'].isin(complex_name)]

comp_regex = re.compile(r'(b[0-9]{4})\((\d?)\)')


def get_coupling_dict(model, mode, atps_name = None, infer_missing_enz=False):
    coupling_dict = get_yeast_coupling_dict(model)
    # coupling_dict.update(get_atp_synthase_coupling(atps_name))
    # if atps_name is not None:
    #     atps = get_atp_synthase_coupling(atps_name)
    #     coupling_dict.update(atps)
    if isinstance(mode, Number):
        kcat = mode
    #elif mode=='kmax':            
    else:
        kcat = get_average_kcat()
    if infer_missing_enz:
        inferred_enz = dict()
        for r in model.reactions:
            if r.id not in coupling_dict \
                    and is_me_compatible(r):
                inferred_enz[r.id] = infer_enzyme_from_gpr(r,
                    default_kcat=get_average_kcat(),
                    default_kdeg=kdeg_enz)

        coupling_dict.update(inferred_enz)
    return coupling_dict


#Velours J, Paumard P, Soubannier V, Spannagel C, Vaillier J, Arselin G,
# Graves PV (May 2000). "Organisation of the yeast ATP synthase F(0):a study
# based on cysteine mutants, thiol modification and cross-linking reagents". 
# Biochimica et Biophysica Acta. 1458 (2–3): 443–56.
# doi:10.1016/S0005-2728(00)00093-1. PMID 10838057
# table 1
atps_general_composition = {
    'YBL099W' :3, # ATP synthase F1 alpha subunit 'ATP1'
    'YJR121W' :3, # ATP synthase F1 beta subunit 'ATP2'
    'YBR039W' :1, # ATP synthase F1 gamma subunit 'ATP3'
    'YDL004W' :1, # ATP synthase F1 delta subunit 'ATP16'
    'YPL271W' :1, # ATP5I  # ATP synthase F1 epsilon subunit 'ATP15'
    'Q0085' :1, # ATP synthase F0 a (6) subunit mitochondrally encoded 'ATP6'
    'Q0080' :1, # ATP synthase F0 A6L (8) subunit mitochondrally encoded 'ATP8' 
    'Q0130' :10, # ATP synthase F0 c subunit mitochondrally encoded 'ATP9'
    'YPL078C' :1, # ATP synthase F0 b subunit 'ATP4'
    'YDR298C' :1, # ATP synthase F0 oscp subunit 'ATP5'
    'YKL016C' :1, # ATP synthase F0 d subunit 'ATP7'
    'YDR377W' :1, # ATP synthase F0 f subunit 'ATP17'
    'YPR020W':1,
    'YML081C-A':1,
    'YLR295C' :1, # ATP synthase F0 h subunit 'ATP14'
} 

def get_atp_synthase_coupling(atps_name):

    # "Between 1404000 per hour at saturation"
    kcat = 1404000*sigma

    def expand_dict(d):
        """
        make one dict per combination of tuple in the keys. Ugly but works.

        :param d:
        :return:
        """
        all_dicts = [dict()]
        for k,v in d.items():
            if isinstance(k,tuple):
                first_choice = k[0]
                new_dicts = []
                for e,this_k in enumerate(k):
                    if e == 0:
                        for this_d in all_dicts:
                            this_d[this_k] = v
                    if e > 0:
                        for this_d in all_dicts:
                            new_d = this_d.copy()
                            new_d[this_k] = new_d.pop(first_choice)
                            new_dicts.append(new_d)
                all_dicts.extend(new_dicts)
            else:
                for this_d in all_dicts:
                    this_d[k] = v

        return all_dicts

    atp_synthases = []

    for e,composition in enumerate(expand_dict(atps_general_composition)):
        this_atp_synthase = Enzyme(atps_name+'_{}'.format(e),
                            name='ATP Synthase_{}'.format(e),
                            kcat_fwd=kcat,
                            kcat_bwd=kcat,
                            kdeg=kdeg_enz,
                            composition=composition)
        atp_synthases.append(this_atp_synthase)

    return {atps_name:atp_synthases}

def get_mrna_dict(model):
    mrna_dict = dict()

    # Generate a mRNA dict

    for x in nt_sequences.index:
        try:
            the_gene = model.genes.get_by_id(x)
        except KeyError:
            model.genes += [cobra.Gene(id=x)]
            the_gene = model.genes.get_by_id(x)

        this_kdeg_mrna = kdeg_mrna # Average value of 5 mins

        new_mrna = mRNA(x,
                        kdeg = this_kdeg_mrna,
                        gene_id = the_gene.id)
        mrna_dict[x] = new_mrna
    return mrna_dict


# Ribosome
#
rrna_genes = ['RDN18-1', 'RDN58-1', 'RDN5-1','RDN25-1']

# Half life of a ribosome is 5 days
# kdeg = ln(2)/t_half * (conversion to hours)
kdeg_rib = np.log(2)/(5*24)

## rPeptides:
# Planta, Rudi J., and Willem H. Mager. "The list of cytoplasmic ribosomal
# proteins of Saccharomyces cerevisiae." Yeast 14.5 (1998): 471-477.
rpeptide_list_a = pd.read_excel(pjoin(data_dir,'cyt_ribosomal_protein_yeast_A.xlsx'),
                             sheet_name='rp',
                             header=0)['Gene'].tolist()
rpeptide_genes_a={k:1 for k in rpeptide_list_a}
rpeptide_genes_a['YBL092W'] = 2

## rPeptides:
# Planta, Rudi J., and Willem H. Mager. "The list of cytoplasmic ribosomal
# proteins of Saccharomyces cerevisiae." Yeast 14.5 (1998): 471-477.
rpeptide_list_b = pd.read_excel(pjoin(data_dir,'cyt_ribosomal_protein_yeast_B.xlsx'),
                             sheet_name='rp',
                             header=0)['Gene'].tolist()
rpeptide_genes_b={k:1 for k in rpeptide_list_b}
rpeptide_genes_b['YBL092W'] = 2

def get_rib():
    """
    # Ribosome

    :return:
    """

    # https://bionumbers.hms.harvard.edu/bionumber.aspx?id=105225&ver=7&trm=Saccharomyces+cerevisiae+translation+rate&org=
    # Bionumber ID 	105225
    # Value         between 10 aa/sec
    
    rib_a = Ribosome(id='rib_a', name='Ribosome_A', kribo=10 * 3600, kdeg=kdeg_rib,
               composition = rpeptide_genes_a, rrna=rrna_genes)
    
    rib_b = Ribosome(id='rib_b', name='Ribosome_B', kribo=10 * 3600, kdeg=kdeg_rib,
               composition = rpeptide_genes_b, rrna=rrna_genes)


    warn('Implement multiple ribosome RNA clusters')
    return rib_a, rib_b

## mitochondrial ribosome
## rPeptides:
# GRAACK, Hanns-Rüdiger, and Brigitte Wittmann-Liebold. 
# "Mitochondrial ribosomal proteins (MRPs) of yeast." Biochemical Journal 329.3 
# (1998): 433-448.
rpeptide_list_mit = pd.read_excel(pjoin(data_dir,'mit_ribosomal_protein_yeast.xlsx'),
                             sheet_name='rp',
                             header=0)['Gene'].tolist()
rpeptide_genes_mit={k:1 for k in rpeptide_list_mit}


#
rrna_genes_mit = ['Q0158', 'Q0020']

# Half life of a ribosome is 5 days
# kdeg = ln(2)/t_half * (conversion to hours)
kdeg_rib = np.log(2)/(5*24)
def get_rib_mit():
    """
    # Ribosome

    :return:
    """

    
    rib = Ribosome(id='rib_mit', name='mit_Ribosome', kribo=10 * 3600, kdeg=kdeg_rib,
               composition = rpeptide_genes_mit, rrna=rrna_genes_mit)


    warn('Implement multiple ribosome RNA clusters')
    return rib

# RNAP
# 12 subunits POLR2A -> POLR2L
rnap_genes = ['YBR154C','YPR187W','YDL140C','YOR151C','YIL021W']
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=109913&ver=6&trm=Saccharomyces+cerevisiae+rna+polymerization+rate&org=
# Bionumber ID  109913
# Value 	    ~ 40 bp/sec

ktrans = 40
def get_rnap():
    
    rnap = RNAPolymerase(id='rnap',
                         name='RNA Polymerase',
                         ktrans = ktrans*3600,
                         kdeg = kdeg_enz,
                         composition = rnap_genes)


    return rnap

# Mitochondrial
rnap_genes_mit = ['YFL036W','YMR228W'] # CPX-3146 in complex portal
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=109913&ver=6&trm=Saccharomyces+cerevisiae+rna+polymerization+rate&org=
# Bionumber ID  109913
# Value 	    ~ 40 bp/sec

# ktrans = 40
def get_rnap_mit():
    
    rnap = RNAPolymerase(id='rnap_mit',
                         name='mitochondrial RNA Polymerase',
                         ktrans = 40*3600,
                         kdeg = kdeg_enz,
                         composition = rnap_genes_mit)


    return rnap

def get_transcription_dict():
    rnap = get_rnap()
    rnap_mit = get_rnap_mit()
    transcription_dict = dict()
    for the_gene_id in gene_list['Gene_id'].values:
        if the_gene_id.startswith('Q'): # mitochondrial gene
            transcription_dict[the_gene_id] = [rnap_mit.id]
        else:
            transcription_dict[the_gene_id] = [rnap.id]
    
    return transcription_dict

def get_translation_dict():
    rib_a, rib_b = get_rib()
    rib_mit = get_rib_mit()
    translation_dict = dict()
    for the_gene_id in pep_list['Gene_id'].values:
        if the_gene_id.startswith('Q'): # mitochondrial gene
            translation_dict[the_gene_id] = [rib_mit.id]
        else:
            translation_dict[the_gene_id] = [rib_a.id, rib_b.id]
    
    return translation_dict

def get_enzyme_fraction(peptides):
    # a function to find the fraction of proteome that is enzyme
    # peptides is the list of model peptides 
    
    # it's possible it cannot find all the enzymes, so a flexibility is introduced
    alpha = 0 # 0% tolerance
    # exclusion = rpeptide_list + rpeptide_list_mit + rnap_genes + rnap_genes_mit
    
    weight_dict_enz = dict()
    for the_peptide in peptides:
        # if the_peptide.id in exclusion:
        #     continue
        weight_dict_enz[the_peptide.id] = [the_peptide.molecular_weight * x for x in \
            prot_abundance[prot_abundance['genes']==the_peptide.id]\
                ['abundance']]
    tot_enz = sum([x[0] for x in weight_dict_enz.values() if len(x)>0 and \
                   not np.isnan(x[0])])
    tot_enz *=1000 #it is in KDa
    
    weight_dict_prot = dict()
    for _, pair in prot_abundance.iterrows():
        gene_id = pair[0]
        abundance = pair[1]
        # if gene_id in exclusion:
        #     continue
        weight_dict_prot[gene_id] = [y * abundance for y in \
                 mws[mws['gene name']==gene_id]['MW']]
    tot_prot = sum([x[0] for x in weight_dict_prot.values() if len(x)>0 and \
                   not np.isnan(x[0])])
    fraction = tot_enz/tot_prot
    return fraction + alpha*fraction

# Setting copy numbers for the genes:
def update_copy_numbers(model):
# Kobayashi, Takehiko. "Regulation of ribosomal RNA gene copy number and its role in modulating genome integrity and evolutionary adaptability in yeast." Cellular and Molecular Life Sciences 68.8 (2011): 1395-1403.
# For rRNAs we need to set the copy number equal to 150. 
    copy_dict = {'RDN18-1' : 150, 'RDN58-1' : 150 ,\
                 'RDN5-1' : 150,'RDN25-1' : 150 ,\
                 'Q0158' : 150, 'Q0020' : 150, }
    for gene,copy_num in copy_dict.items():
        try:
            this_gene = model.genes.get_by_id(gene)
        except KeyError:
            continue
        this_gene.copy_number = copy_num
    return

def trna_charging_enz(model):
# a function to define enzymes for tRNA charging reactions    
    aa_dict, _, _, _ = get_monomers_dict()
    
    annot_dict = {'A' : 'YOR335C',
                  'R' : 'YDR341C',
                  'N' : 'YHR019C',
                  'D' : 'YLL018C',
                  'C' : 'YNL247W',
                  'E' : 'YOR168W',
                  'Q' : 'YGL245W',
                  'G' : '( YBR121C or YPR081C )',
                  'I' : 'YBL076C',
                  'L' : 'YPL160W',
                  'K' : 'YDR037W',
                  'M' : 'YGR264C',
                  'F' : '( YFL022C and YLR060W )',
                  'P' : 'YHR020W',
                  'S' : '( YDR023W or YHR011W )',
                  'T' : 'YIL078W',
                  'W' : 'YOL097C',
                  'Y' : 'YGR185C',
                  'H' : 'YPR033C',
                  'V' : 'YGR094W'}
    
    enzyme_list = list() # a list for adding all enzymes to the model
    enzyme_dict = dict() # a dict to keep the relation between enzymes and rxns
    for aa, gpr in annot_dict.items():
        compositions = compositions_from_string(model, gpr)
        new_enzymes = list()
        for e,composition in enumerate(compositions):
            new_enzyme = Enzyme(id = 'trna_ch_{}_{}'.format(aa_dict[aa],e),
                                kcat=Kcat_average,
                                kdeg=kdeg_enz,
                                composition=composition)
            new_enzymes.append(new_enzyme)
        enzyme_list += new_enzymes
        enzyme_dict[aa_dict[aa]] = new_enzymes
    return enzyme_list, enzyme_dict

def remove_DNA(model):
    # if allocation data is available we can remove DNA from biomass, at it 
    # is made growth dependent in the code.
    DNA = model.metabolites.get_by_id('s_3720_c')
    DNA_rxn = model.reactions.get_by_id('r_4050')
    model.remove_metabolites(DNA)
    model.remove_reactions(DNA_rxn)

# for some procedures we need to have unmodified model
cobra_model = cobra.io.load_matlab_model('../models/yeast8_thermo_curated.mat')
    
def get_macromole_ratio():
    # model must be an unmodified model
    mass_ratios = find_bbb_ratio(cobra_model)
    
    # mass_ratios = {'protein': 0.4448416049408599,
    #  'RNA': 0.07292485326899344,
    #  'DNA': 0.00416713447251391,
    #  'lipid': 0.06667415156022256,
    #  'carbohydrate': 0.422964148960162,
    #  'ion': 0.029169941307597376,
    #  'cofactor': 0.0052089180906423884,
    #  'total mass': 1.0459507526009915}
    return mass_ratios

def modify_GAM(model, growth_id, prot_rxn_id=None):
    '''
    This function modify GAM by removing the energetic cost of peptide synthesis.

    Parameters
    ----------
    model: cobra-model, with modified biomass
    growth_id: str
    prot_rxn_id : str, optional
        The id of protein pseudoreaction if aa are not taken into account in the
        biomass reaction. The default is None.

    Returns
    -------
    cobra-model with modified GAM

    '''
    # model must be an unmodified model
    growth_rxn = cobra_model.reactions.get_by_id(growth_id)
    
    if prot_rxn_id is None:
         aa_dict,_,_,_ = get_monomers_dict()
         aa_stoic = [abs(growth_rxn.get_coefficient(x)) \
                     for _,x in aa_dict.items()]  
    else:
         prot_rxn = cobra_model.reactions.get_by_id(prot_rxn_id)
         aa = [x.id for x in prot_rxn.reactants]
         aa_stoic = [abs(prot_rxn.get_coefficient(x)) \
                     for x in aa]
     
    tot_aa = sum(aa_stoic) # number of moles required to synthesize one unit of biomass
    gtp_expense = 2 * tot_aa # 2 moles of GTP needed per each mole of aa
    # now we should find GAM metabolites
    GAM_mets = {'atp':-1, 'h2o':-1,
                'adp':1,'h':1,'pi':1} # differentiating reactants & products
    essentials = get_essentials()
    # create a dict to subtract from metabolite ids
    mets = {essentials[k]:v*gtp_expense for k,v in GAM_mets.items()}
    # to apply changes on the modified model
    biomass = model.reactions.get_by_id(growth_id)
    biomass.subtract_metabolites(mets)
    
def constrain_enzymes(model, mass_ratios):
    # a function to add a constraint on total amount of enzymes based on their
    # fraction from total amount of proteins (should be before adding dummy)
    enz_ratio = get_enzyme_fraction(model.peptides) 
                
    enz_vars = model.get_variables_of_type(EnzymeVariable)
    
    # we should first exclude dummy, ribosomes and rnaps
    exclusion = ['dummy_enzyme', #'rib', 'rib_mit', 'rnap', 'rnap_mit'
                 ]
    exclusion = ['EZ_{}'.format(x) for x in exclusion]
    enz_vars = [x for x in enz_vars if x.name not in exclusion]
    
    expr = symbol_sum([x for x in enz_vars])
    try:
        # modle with variable biomass
        model.add_constraint(kind = ConstantAllocation, 
                                 hook = model, 
                                 expr = expr - enz_ratio * model.variables.IV_prot_ggdw,
                                 id_ = 'enzyme_fix',
                                 ub = 0)
    except AttributeError:
        model.add_constraint(kind = ConstantAllocation, 
                                 hook = model, 
                                 expr = expr,
                                 id_ = 'enzyme_fix',
                                 lb = 0, # cannot be negative
                                 ub = enz_ratio * mass_ratios['protein'])
