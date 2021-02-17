from copy import deepcopy
import numpy as np

from ..optim.variables import  BinaryVariable, \
                                     GeneVariable, get_binary_type
from ..optim.constraints import RNAPAllocation

import sympy
import pandas as pd

from ..utils.parsing import parse_gpr


def is_transport(rxn):
    # taking a reaction ID, it determines if it's a trasport or not
    # transport reaction is a reaction that just move a metabolite from one compartmenst to another
    rmets = rxn.metabolites.keys()
    if len(rmets) == 2:
        comps = [met.compartment for met in rmets]
        if len(comps) > 1:
            return True
        else:
            return False
    else:
        return False   

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


class RNAPActivator(GeneVariable,BinaryVariable):
    """
    Class to represent a binary variable that activates a catalytic constraint
    or relaxes it
    """
    def __init__(self, gene, **kwargs):
        GeneVariable.__init__(self, gene=gene,
                                  type = get_binary_type(),
                                  **kwargs)

    prefix = 'RPA_'

def relax_rnap_constraints(model, min_growth):
    """
    Find a minimal set of rnap_allocation constraints to relax to meet a minimum
    growth criterion

    :param model:
    :type model: :class:`etfl.core.memodel.MEModel`
    :param min_growth:
    :return:
    """

    relaxation = deepcopy(model)

    relaxation.objective = relaxation.growth_reaction

    new_objective = 0

    for the_constr in relaxation.rnap_allocation:
        activator = relaxation.add_variable(kind = RNAPActivator,
                                               hook = the_constr.gene,
                                               )
        
        new_expr = the_constr.constraint.expression - (1-activator) * relaxation.big_M

        lb = the_constr.constraint.lb
        ub = the_constr.constraint.ub

        relaxation.remove_constraint(the_constr)
        relaxation.add_constraint(kind = RNAPAllocation,
                                     hook = the_constr.gene,
                                     expr = new_expr,
                                     lb = lb,
                                     ub = ub)
        


        new_objective += activator

    relaxation.repair()

    relaxation.growth_reaction.lower_bound = min_growth

    relaxation.objective = new_objective

    relaxation.optimize()

    activators = relaxation.get_variables_of_type(RNAPActivator)
    activator_states=relaxation.solution.raw.loc[activators.list_attr('name')]

    for the_var in activators.list_attr("variable"):
        the_var.lb = int(relaxation.solution.raw.loc[the_var.name])
        the_var.ub = int(relaxation.solution.raw.loc[the_var.name])

    relaxed_model = deepcopy(model)

    for act in activators:
        if np.isclose(activator_states[act.name],0):
            the_cons = relaxed_model.forward_catalytic_constraint.get_by_id(act.id)
            relaxed_model.remove_constraint(the_cons)

    return activator_states, relaxed_model, relaxation


def relax_alloc_cstr(model):
    # a function to relax rnap_allocation constraints for rnap & rib
    cstr_set = model.get_constraints_of_type(RNAPAllocation)
    relax_set = [x for x in cstr_set \
                 if x.id in model.ribosome['rib'].composition or \
                 x.id in model.rnap['rnap'].composition]
    for x in relax_set:
        model.remove_constraint(x)
        
def relax_thermo(model):
        # a function to relax thermo constraints for some reactions
        rxn_ids = ['r_1758','r_3552','r_3721','r_3722'] # for producin ergosterol (lipidplots)
        rxns = [model.reactions.get_by_id(x) for x in rxn_ids]
        for rxn in rxns:
            rxn.thermo['computed'] = False
        
def enzrxn_to_gpr(rxn):
    try:
        z = rxn.enzymes
    except (AttributeError):
            return 
    return ' OR '.join((' AND '.join(['{}*{}'.format(v,k) for v,k \
                                      in isozyme.composition.items()])) \
                                        for isozyme in rxn.enzymes)

def compositions_from_string(model, this_gpr):
    """
    *Warning*: Use this function only if you have no information on the enzymes.
    Logically parses the GPR to automatically find isozymes ( logical OR )
    and subunits ( logical AND ), and creates the necessary complexation
    reactions: 1 per isozyme, requiring the peptides of each subunit

    :param reaction:
    :type reaction: cobra.Reaction
    :return:
    """


    sym_gpr = parse_gpr(this_gpr)

    if isinstance(sym_gpr, sympy.Symbol):
        # GPR of the type: '(gene0)'
        # Gene <=> Protein
        isozymes = [sym_gpr]
    elif isinstance(sym_gpr, sympy.And):
        # GPR of the type: '(gene0 & gene1)'
        # Subunits of one enzyme
        isozymes = [sym_gpr]
    elif isinstance(sym_gpr, sympy.Or):
        # GPR of the type: '(gene0 | gene1)', '((gene0 & gene1) | gene2)'
        # Two isozymes that are the arguments of the OR
        isozymes = sym_gpr.args

    compositions = []

    for e, this_isozyme in enumerate(isozymes):
        if isinstance(this_isozyme, sympy.And):
            # this is a GPR with several subunits
            peptides = {x.name.upper().replace('_','-'): 1 \
                        for x in this_isozyme.args}
        elif isinstance(this_isozyme, sympy.Symbol):
            # there is only one subunit
            peptides = {this_isozyme.name.upper().replace('_','-'): 1}
        else:
            # The GPR has been incorrectly parsed
            model.logger.error('Incorrect parsing of {}'.format(isozymes))
            raise TypeError

        compositions += [peptides]

    return compositions

def special_coupling(rxn, enz_dict, genes_dict):
    # a function to handle specific cases for enzyme-reaction coupling,
    # when at least a gene is common between isozymes
    
    # Let's first split the isozymes
    gpr = rxn.gene_reaction_rule
    enz_list = gpr.replace('(','').replace(')','').split(' or ')
    matched_enz = []
    for this_enz in enz_list:
        gene_ids = this_enz.strip().split(' and ')
        score_df = pd.Series({k:jaccard_score(gene_ids,v) for k,v in genes_dict.items()})
        perfect_scores = score_df[score_df > 0.91]
        matched_enz += [enz_dict[x] for x in perfect_scores.index]
    
    return matched_enz

def coupling_trna_enzymes_dict(model, enz_dict):
    couple_dict = dict()
    rxn_id_mod = 'trna_ch_'
    for aa, enzymes in enz_dict.items():
        this_rxn = model.reactions.get_by_id(rxn_id_mod + aa)
        couple_dict[this_rxn.id] = enzymes
        
    return couple_dict

def removing_excessive_mets_rxns(model):
    
    trna_charged_mets = ['s_0404_c',
    's_0428_c',
    's_0430_c',
    's_0432_c',
    's_0542_c',
    's_0747_c',
    's_0748_c',
    's_0757_c',
    's_0832_c',
    's_0847_c',
    's_1077_c',
    's_1099_c',
    's_1148_c',
    's_1314_c',
    's_1379_c',
    's_1428_c',
    's_1491_c',
    's_1527_c',
    's_1533_c',
    's_1561_c',]
    
    trna_charged_rxns = ['r_0157',
    'r_0209',
    'r_0212',
    'r_0220',
    'r_0313',
    'r_0478',
    'r_0479',
    'r_0512',
    'r_0539',
    'r_0665',
    'r_0701',
    'r_0711',
    'r_0729',
    'r_0852',
    'r_0941',
    'r_0995',
    'r_1042',
    'r_1057',
    'r_1066',
    'r_1089',]
    
    pseudo_rxns = ['r_4047','r_4049',] #'r_4050']
    pseudo_mets = ['s_3717_c','s_3719_c',] #'s_3720_c']
    
    rxn_list = pseudo_rxns + trna_charged_rxns
    met_list = pseudo_mets + trna_charged_mets
    
    rxn_set = [model.reactions.get_by_id(x) for x in rxn_list]
    met_set = [model.metabolites.get_by_id(x) for x in met_list]
    
    for rxn in rxn_set:
        model.remove_reactions(rxn)
    for met in met_set:
        model.remove_metabolites(met)
    return
        
def find_bbb_ratio(model, BBB = ['all']):
    '''
    This function is designed for Yeast8!
    
    A function to find the mass ration of each biomass building block in FBA biomass definition.
    Also can return total weight of biomass
    inputs:
       model: a cobra model before any modification
       BBB: A list of the building block(s) of interest [lipid, carbohydrate, protein, RNA, DNA, ion, cofactor, {all}]
    '''
    
    growth_rxn = 'r_4041'
    pseudo_rxns = {'protein':'r_4047',
    'carbohydrate':'r_4048',
    'RNA':'r_4049',
    'DNA':'r_4050',
    'lipid backbone':'r_4063',
    'lipid chain':'r_4065',
    'cofactor':'r_4598',            
    'ion':'r_4599',}
    
    pseudo_mets = {'s_1096_c':'lipid', 
                     's_3717_c':'protein',
                     's_3718_c':'carbohydrate',
                     's_3719_c':'RNA',  
                     's_3720_c':'DNA',
                     's_4205_c':'cofactor',
                     's_4206_c':'ion',}
    
    ratios = dict()
    
    if 'all' in BBB:
        BBB = ['lipid', 'carbohydrate', 'protein', 'RNA', 'DNA', 'ion', 'cofactor','all']
    if 'protein' in BBB:
        pseudo_id = pseudo_rxns['protein']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        # the protein psedoreaction is junk. I should find AAs first.
        AA_met = {'A': 's_0955_c',
               'R': 's_0965_c',
               'N': 's_0969_c',
               'D': 's_0973_c',
               'C': 's_0981_c',
               'E': 's_0991_c',
               'Q': 's_0999_c',
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
               'W': 's_1048_c',
               'Y': 's_1051_c',
               'V': 's_1056_c', }
        AA_MWs = {k:model.metabolites.get_by_id(v).formula_weight \
                  for k,v in AA_met.items()}
        tRNA_met = {'A' : 's_0404_c',
                    'R' : 's_0428_c',
                    'N' : 's_0430_c',
                    'D' : 's_0432_c',
                    'C' : 's_0542_c',
                    'E' : 's_0747_c',
                    'Q' : 's_0748_c',
                    'G' : 's_0757_c',
                    'H' : 's_0832_c',
                    'I' : 's_0847_c',
                    'L' : 's_1077_c',
                    'K' : 's_1099_c',
                    'M' : 's_1148_c',
                    'F' : 's_1314_c',
                    'P' : 's_1379_c',
                    'S' : 's_1428_c',
                    'T' : 's_1491_c',
                    'W' : 's_1527_c',
                    'Y' : 's_1533_c',
                    'V' : 's_1561_c'}
        AA_stoic = {k:pseudo_rxn.get_coefficient(v) for k,v in tRNA_met.items()}
        prot_ratio = sum([-v * AA_MWs[k] for k,v in AA_stoic.items()])/1000
        ratios['protein'] = prot_ratio
    if 'RNA' in BBB:
        pseudo_id = pseudo_rxns['RNA']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        NMP_met = {
                'u': 's_1545_c',
                'g': 's_0782_c',
                'a': 's_0423_c',
                'c': 's_0526_c'}
        NMP_MWs = {k:model.metabolites.get_by_id(v).formula_weight \
                  for k,v in NMP_met.items()}
        NMP_stoic = {k:pseudo_rxn.get_coefficient(v) for k,v in NMP_met.items()}
        RNA_ratio = sum([-v * NMP_MWs[k] for k,v in NMP_stoic.items()])/1000
        ratios['RNA'] = RNA_ratio
    if 'DNA' in BBB:
        pseudo_id = pseudo_rxns['DNA']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        dNMP_met = {
                't': 's_0649_c',
                'g': 's_0615_c',
                'a': 's_0584_c',
                'c': 's_0589_c'}
        dNMP_MWs = {k:model.metabolites.get_by_id(v).formula_weight \
                  for k,v in dNMP_met.items()}
        dNMP_stoic = {k:pseudo_rxn.get_coefficient(v) for k,v in dNMP_met.items()}
        DNA_ratio = sum([-v * dNMP_MWs[k] for k,v in dNMP_stoic.items()])/1000
        ratios['DNA'] = DNA_ratio
    if 'lipid' in BBB:
        pseudo_id = pseudo_rxns['lipid chain']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        LC_MWs = {met.name:met.formula_weight for met in pseudo_rxn.metabolites.keys()}
        LC_stoic = {met.name:pseudo_rxn.get_coefficient(met.id) \
                       for met in pseudo_rxn.metabolites.keys()}
        
        pseudo_id = pseudo_rxns['lipid backbone']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        LB_MWs = {met.name:met.formula_weight for met in pseudo_rxn.metabolites.keys()}
        # fatty acid backbone doesn't have MW; its stoichiometric coeff is 0.0014801
        # Based on producing reactions its averaged MW is 730
        LB_MWs['fatty acid backbone [cytoplasm]'] = 730
        
        LB_stoic = {met.name:pseudo_rxn.get_coefficient(met.id) \
                       for met in pseudo_rxn.metabolites.keys()}
        
        lipid_ratio = sum([-v * LC_MWs[k] for k,v in LC_stoic.items() if v<0])/1000 + \
                        sum([-v * LB_MWs[k] for k,v in LB_stoic.items() if v<0])/1000
        ratios['lipid'] = lipid_ratio
    if 'carbohydrate' in BBB:
        pseudo_id = pseudo_rxns['carbohydrate']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        sugar_MWs = {met.name:met.formula_weight for met in pseudo_rxn.metabolites.keys()}
        sugar_stoic = {met.name:pseudo_rxn.get_coefficient(met.id) \
                       for met in pseudo_rxn.metabolites.keys()}
        carbohydrate_ratio = sum([-v * sugar_MWs[k] for k,v in sugar_stoic.items() if v<0])/1000
        ratios['carbohydrate'] = carbohydrate_ratio
    if 'ion' in BBB:
        pseudo_id = pseudo_rxns['ion']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        ions_MWs = {met.name:met.formula_weight for met in pseudo_rxn.metabolites.keys()}
        ions_stoic = {met.name:pseudo_rxn.get_coefficient(met.id) \
                       for met in pseudo_rxn.metabolites.keys()}
        ion_ratio = sum([-v * ions_MWs[k] for k,v in ions_stoic.items() if v<0])/1000
        ratios['ion'] = ion_ratio
    if 'cofactor' in BBB:
        pseudo_id = pseudo_rxns['cofactor']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        cofactor_MWs = {met.name:met.formula_weight for met in pseudo_rxn.metabolites.keys()}
        cofactor_stoic = {met.name:pseudo_rxn.get_coefficient(met.id) \
                       for met in pseudo_rxn.metabolites.keys()}
        cofactor_ratio = sum([-v * cofactor_MWs[k] for k,v in cofactor_stoic.items() if v<0])/1000
        ratios['cofactor'] = cofactor_ratio
    
    pseudo_rxn = model.reactions.get_by_id(growth_rxn)
    for k,v in pseudo_mets.items():
        ratios[v] *= abs(pseudo_rxn.get_coefficient(k))
    if 'all' in BBB:
        pseudo_rxn = model.reactions.get_by_id(growth_rxn)
        ratios['total mass'] = sum([x for x in ratios.values()])
    
    return ratios

def stoich_coeff_from_mass_ratio(model, alloc_data, macromole_ratio, \
                                 growth_id, biomass_comps, met_macromol_dict, \
                                 lumped_comp = True):
    '''
    This is a function to find new stoichiometric coefficients for biomass reaction
    based on experimental data for mass ratios.
    
    inputs:
        model: a cobra or pytfa model (not modified)
        alloc_data: a dict or DataFrame of different biomass components with
            their experimental mass ratios
        macromole_ratio: dict of mass ratios of different biomass components in
            GEM.
        growth_id: the rxn id for growth.
        biomass_comps: a list of biomass components (metabolite ids),
        lumped_comp: indicates if biomass metabolites are lumped into pseudometabolites
        met_macromol_dict: adictionary to relate metabolite ids with each type
            of macromolecules, e.g. protein, RNA, etc. Keys are met ids and values
            are macromolecule names compatible with tags in alloc_data & macromole_ratio
        
    ouputs:
        new model with modified stoichiometric coefficients in its growth reaction.
    '''
    
    # finding the original stoichiometric coefficients for each metabolite
    biomass_rxn = model.reactions.get_by_id(growth_id)
    
    if lumped_comp:
        # nothing to do :)
        org_coeffs = {x : biomass_rxn.get_coefficient(x) for x in biomass_comps}
    else:
        # first, we should lump different metabolites to relate them to macromolecules
        raise NotImplementedError()
        
    try:
        tot_mass = macromole_ratio['total mass']
    except KeyError:
        tot_mass = sum([x for x in macromole_ratio.values()])
        
    ratios = {k : v/tot_mass for k,v in macromole_ratio.items()}
    new_coeffs = {k : v * alloc_data[met_macromol_dict[k]] /\
                  ratios[met_macromol_dict[k]] for k,v in org_coeffs.items()}
    
    # the new stoichiometric coefficient is added to the previous one
    # to avoid redundancy I should first subtract the old coefficient
    change_coeffs = {k : v - org_coeffs[k] for k,v in new_coeffs.items()}
    biomass_rxn.add_metabolites(change_coeffs)
    
    return  

def _infer_manual_information(df):
    kcat_dict = dict()
    genes_dict = dict()
    rxn_dict = dict()

    for ind ,row in df.iterrows():
        this_id = 'enzyme_{}'.format(ind)
        gene_name = row['composition']
        kcat = row['kcat']
        rxn = row['rxnId']
        rxn_dict[this_id] = rxn
        kcat_dict[this_id] = kcat
        genes_dict[this_id]= tuple(gene_name.split('//'))
    return genes_dict, rxn_dict, kcat_dict

def extract_kcat_from_gecko(model): 
    # model is needed just to extract enzyme compositions 
    # I need a function to infer kcat data from gecko.
    # the output of matlab script is a list of enzyme(composition), rxn_id, kcat
    # should be converted to a list of enzyme_rxn_id, kcat
    kcat_automatic          = pd.read_excel('../data/enzyme2kcat.xlsx',
                                      sheet_name = 'enzyme2kcat')
    kcat_manual          = pd.read_excel('../data/gecko_enzyme_rxn_pairs.xls',
                                      sheet_name = 'Sheet1')
    kcat_manual['kcat'] *= 3600
    found = [] # list oof indices of gecko enzymes that are found in yETFL
    genes_dict, rxn_dict, kcat_dict = _infer_manual_information(kcat_manual)
    
    final_dict = dict()
    for enz in model.enzymes:
        genes_of_enzyme = [g for g in enz.composition.keys()]
        score_df = pd.Series({k:jaccard_score(genes_of_enzyme,v) \
                              for k,v in genes_dict.items()})
        # top_scores = score_df.max()
        perfect_scores = score_df[ score_df > 0.4]
        if len(perfect_scores) > 0: # has a match in manual kcats
            for ind in perfect_scores.index:
                if not isinstance(rxn_dict[ind],str): # usual value for kcat
                    final_dict[enz.id] = kcat_dict[ind]
                else: # exceptional value for kcat 
                    final_dict[enz.id + '_' + rxn_dict[ind]] = kcat_dict[ind]
                found.append(ind)
        else: # without any match in manual kcats
            kcat = kcat_automatic[kcat_automatic['complex'] == enz.id]['kcat']
            final_dict[enz.id] = float(kcat)
    
    table = pd.DataFrame(final_dict.items(), columns=['complex','kcat'])
    table.to_csv('../data/complex2kcat.csv')
    return final_dict, found   

def infer_composition(unfound_enzymes):
    # a function to create composition database for unfound enzymes by assuming a stoichiometric coefficient of 1
    # unfound enzymes is a modified output of extract_kcat_from_gecko(model)
    gene_sets = [str(x) for x in unfound_enzymes['0']]
    kcats = [str(x) for x in unfound_enzymes['1']] # to create kcat list with similar IDs
    # we need to create four lists corresponding to each column of compsoition data
    ids = []
    gene_ids = []
    stoichiometry = []
    products = []
    kcat_dict = dict()
    for idx, genes in enumerate(gene_sets):
        this_id = 'COMPLEX_{}'.format(idx)
        gene_names = genes.split('//')
        bbb1 = gene_names
        bbb2 = ['TUPLE // {}-MONOMER // 1'.format(x) for x in gene_names]
        bbb3 = ['{}-MONOMER'.format(x) for x in gene_names]
        ids += [this_id]
        gene_ids += [' // '.join(bbb1)]
        stoichiometry += [' // '.join(bbb2)]
        products += [' // '.join(bbb3)]
        kcat_dict[this_id] = kcats[idx]
    
    compositions = pd.DataFrame([ids, gene_ids, stoichiometry, products], \
        index = ['Product Name','Genes','Component coefficients','Gene products'])
    kcat_data = pd.DataFrame(kcat_dict.items())
    return compositions, kcat_data
        
        
    
    
    

