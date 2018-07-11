import cobra
import pandas  as pd
import re

from collections import defaultdict

B_GENE_REGEX = re.compile(r'b[0-9]{4}')
EC_NUMBER_REGEX_HINTS = re.compile(r'EC [0-9]+\.[0-9]+\.[0-9]+\.[0-9]+')
SABIORK_REGEX = re.compile(r'sabiork:[0-9]+')

info_peptides  = pd.read_csv('info_ecoli/bigg_iJO1366.peptides.tsv',
                             header = 0,
                             delimiter = '\t')
info_peptides.columns = ['peptide','description','xrefs','gene']
info_reactions = pd.read_csv('info_ecoli/bigg_iJO1366.reactions.tsv',
                             header = 0,
                             delimiter = '\t')
info_reactions.columns = ['sid','reac','source','equa','EC','pathway','xref']

def find_b_gene(gene_list):
    try:
        return B_GENE_REGEX.findall(gene_list)[0]
    except IndexError:
        return None

def find_ec_hints(desc):
    return EC_NUMBER_REGEX_HINTS.findall(desc)

def listify_ecs(ec_str):
    return str(ec_str).split(';')

def find_sabiork(exref_str):
    return SABIORK_REGEX.findall(str(exref_str))

info_peptides['b_gene'] = info_peptides['gene'].apply(find_b_gene)
info_peptides['ec_hints'] = info_peptides['description'].apply(find_ec_hints)

info_reactions['ec_list'] = info_reactions['EC'].apply(listify_ecs)
info_reactions['sabiork'] = info_reactions['xref'].apply(find_sabiork)

ec_cobra = cobra.io.load_json_model('../../models/iJO1366.json')# Make it human_readable

def match_cobra_rxn_to_mnx(rid):
    return info_reactions[info_reactions['source'] == 'R_'+rid]

for reaction in ec_cobra.reactions:
    row = match_cobra_rxn_to_mnx(reaction.id)
    try:
        reaction.notes['ec_numbers'] = row['ec_list'].tolist()[0]
        reaction.notes['sabiork'] = row['sabiork'].tolist()[0]
    except IndexError:
        pass

for gene in ec_cobra.genes:
    try:
        gene.notes['ec_hints'] = info_peptides[info_peptides['b_gene'] == gene.id]['ec_hints'].tolist()[0]
    except IndexError:
        pass

cobra.io.save_json_model(ec_cobra, 'iJO1366_with_xrefs.json')