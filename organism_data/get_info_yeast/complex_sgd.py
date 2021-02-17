#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Inferring protein complex information from SGD database
Complex Portal: https://www.ebi.ac.uk/complexportal/complex/organisms
"""
import xml.etree.ElementTree as ET
import pandas as pd
from warnings import warn


GENE_IDS=[]
COMP_COEF=[]
PROD_IDS=[]
COMPLEX_NAME=[]
for x in range(0, 580):
    data_file='yeast_complex_{}.xml'.format(x+1)
    tree=ET.parse(data_file)
    iterator=tree.getiterator()
    comm_sys=dict() # converting common gene names to systematic gene names
    parent_map = dict((c, p) for p in iterator for c in p)  # mapping each node to its parents
     
    
    for i in iterator:
        if ('type','complex systematic name') in i.items():
            Stoich=i.text.__str__()
            Stoich=Stoich.replace(' ','') # deleting whitespaces
        elif ('type','locus name') in i.items():
            parent=parent_map[i]
            for j in list(parent):   # iterating in sibling nodes
                if ('type','gene name') in j.items():
                    comm_sys[j.text.__str__()]=i.text.__str__()
        elif ('db','complex portal') in i.items():
            complex_name=i.get('id')
            
    # Sometimes 2x[GEN1,GEN2] structure is used
    if '[' in Stoich:
        if Stoich[0].isdigit():
            helper=Stoich.split('[')
            gene_set=helper[1].split(']')[0]
            gene_set=  gene_set.split(':')
            composition = ['{}{}'.format(helper[0],g) for g in gene_set]
        else:
            warn('an unknown structure encountered')
            print(x)
    elif ':' in Stoich:
        composition = Stoich.split(':')
    elif ',' in Stoich:
        composition = Stoich.split(',') # genes are separeted sometimes by comma somtimes by colon
    
    if len(composition) - len(comm_sys) != 0:   # This means the composition is not well_known
        continue
    
    gene_ids=[]
    comp_coef=[]
    prod_ids=[]
    for i in composition:
        coef=i.split('x')
        # Sometimes, they've used X for splitting!
        if len(coef)==1:
            for s in coef:
                if s[0].isdigit() and s[1]=='X':
                    coef=i.split('X')
        # coef[0] is stoichiometry coefficient and coef[1] is common gene name
        if len(coef) == 1: # we should define a different rule for components with 1 coefficient in stoichiometry
            coef += coef
            coef[0]=1
        coef[1]=comm_sys[coef[1].upper()]
        gene_ids=gene_ids + [coef[1]]
        comp_coef=comp_coef + ['TUPLE // {}-MONOMER // {}'.format(coef[1],coef[0])]
        prod_ids=prod_ids + ['{}-MONOMER'.format(coef[1])]
        
    gene_ids = ' // '.join(gene_ids)
    comp_coef = ' // '.join(comp_coef)
    prod_ids = ' // '.join(prod_ids)
    GENE_IDS = GENE_IDS + [gene_ids]
    COMP_COEF = COMP_COEF + [comp_coef]
    PROD_IDS = PROD_IDS + [prod_ids]
    COMPLEX_NAME = COMPLEX_NAME + [complex_name]

Final_result=pd.DataFrame({'Product Name': COMPLEX_NAME, 'Genes': GENE_IDS, 
                           'Component coefficients': COMP_COEF, 'Gene products': PROD_IDS})    
Final_result.to_csv('Yeast_complex_portal.csv', )    

