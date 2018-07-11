import pandas as pd
import numpy as np
import re

raw_cpx2pep  = pd.read_csv('info_ecoli/ecocyc_complex2peptides.csv')
raw_cpx2pep.columns = ['complex','putative_genes','coeffs','peptides','var4','var5','var6']

raw_pep2gene = pd.read_csv('info_ecoli/ecocyc_peptides2genes.csv')
raw_pep2gene.columns = ['peptide','monomer_id','gene']

raw_cpx2ec   = pd.read_csv('info_ecoli/ecocyc_complex2EC.csv')
raw_cpx2ec = raw_cpx2ec[raw_cpx2ec.columns[:2]]
raw_cpx2ec.columns = ['complex', 'ec']

ECO_TUPLE = re.compile(r"(TUPLE \/\/)\s([A-Z0-9\-]*)\s(\/\/)\s(\d)",
                       flags = re.I)


def listify(df, lst_cols, delimiter = ','):
    for col in lst_cols:
        df[col] = df[col].str.split(delimiter)
    return df

def remake_tuples(row):
    all_eco_tuples = ECO_TUPLE.findall(row)
    return pd.Series([[int(x[3]) for x in all_eco_tuples], # return coefficient
                      [    x[1]  for x in all_eco_tuples]]) # return monomer id
    # return [{x[1]:x[3]} for x in all_eco_tuples]


def find_inconsistent(df, cols):
    n = len(cols)
    li = [0]*n
    for e,col in enumerate(cols):
        li[e] = df[col].apply(len)

    bool_table = pd.concat( [li[x] == li[y] for x in range(n) for y in range(x)],
                            axis=1)
    consistent   = df[ bool_table.apply(all,axis=1)]
    inconsistent = df[~bool_table.apply(all,axis=1)]
    return consistent, inconsistent

def explode(df, lst_cols, fill_value=''):
    """
    Inspired from
    https://stackoverflow.com/questions/12680754/split-explode-pandas-dataframe-string-entry-to-separate-rows

    :param df:
    :return:
    """

    # make sure `lst_cols` is a list
    if lst_cols and not isinstance(lst_cols, list):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)

    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()

    if (lens > 0).all():
        # ALL lists in cells aren't empty
        return pd.DataFrame(
            {col: np.repeat(df[col].values, df[lst_cols[0]].str.len()) for col
             in idx_cols}).assign(
            **{col: np.concatenate(df[col].values) for col in lst_cols}).loc[:,
               df.columns]
    else:
        # at least one list in cells is empty
        return pd.DataFrame(
            {col: np.repeat(df[col].values, df[lst_cols[0]].str.len()) for col
             in idx_cols}).assign(
            **{col: np.concatenate(df[col].values) for col in lst_cols}).append(
            df.loc[lens == 0, idx_cols]).fillna
        (fill_value).loc[:, df.columns]


DELIMITER = ' // '


cols_cpx2pep = ['peptides','coeffs']

tmp_cpx2pep = raw_cpx2pep.dropna(subset=cols_cpx2pep)
tmp_cpx2pep[['coeffs','obj_ids']] = tmp_cpx2pep['coeffs'].apply(remake_tuples)
tmp_cpx2pep = listify(tmp_cpx2pep,
                       lst_cols=['peptides'],
                       delimiter=DELIMITER)

tmp_cpx2pep, inconsistent_cpx2pep = find_inconsistent(tmp_cpx2pep,
                                                      cols=cols_cpx2pep) # genes do not coincide with subunits

cpx2pep     = explode(tmp_cpx2pep.reset_index(),
                             lst_cols=cols_cpx2pep+['obj_ids'],
                             )

tmp_cpx2ec = raw_cpx2ec.dropna(subset=['ec'])
tmp_cpx2ec = listify(tmp_cpx2ec,
                     lst_cols=['ec'],
                     delimiter=DELIMITER)
cpx2ec = explode(tmp_cpx2ec.reset_index(),
                 lst_cols=['ec'])
cpx2ec['ec'] = cpx2ec['ec'].str.replace('EC-','')

# cpx2pep[['complex','genes','coeffs','obj_ids']].to_csv('info_ecoli/cleaned_up_cpx2pep.csv')

df1 = cpx2pep[['complex','putative_genes','coeffs','obj_ids']]
result = df1.merge(raw_pep2gene, how='left', left_on='obj_ids', right_on='monomer_id')

result[['complex',
        'putative_genes',
        'coeffs','monomer_id',
        'gene','obj_ids']].to_csv('info_ecoli/complex2genes.csv')

cpx2ec[['complex','ec']].to_csv('info_ecoli/complex2ec.csv')



