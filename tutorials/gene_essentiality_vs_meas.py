import pandas as pd
import numpy as np
from numpy import linspace
from sklearn.metrics import matthews_corrcoef, confusion_matrix

read_ijo_essentiality_data = lambda x:pd.read_excel(
                                    '../organism_data/info_ecoli/orth2011_ijo1366_SI.xlsx',
                                    sheet_name = 'Table 13',
                                    # header = 2,
                                    usecols = x).dropna()

ijo_tp = read_ijo_essentiality_data('O').iloc[2:] #drop extra lines
ijo_fp = read_ijo_essentiality_data('P').iloc[1:]
ijo_fn = read_ijo_essentiality_data('Q').iloc[1:]
ijo_tn = read_ijo_essentiality_data('R').iloc[1:]

for df in [ijo_tp, ijo_fp, ijo_fn, ijo_tn]:
    df.columns = ['genes']

# According to Orth et al 2011
# Non essential genes are TP + FN
non_ess_genes = pd.concat([ijo_tp,ijo_fn], axis=0, ignore_index = True)

# Essential genes are TN + FP
essential_genes = pd.concat([ijo_tn, ijo_fp], axis=0, ignore_index=True)

# etfl_pred = pd.read_csv('outputs/gene_essentiality_vETFL65.csv', index_col = 0)
# etfl_pred = pd.read_csv('outputs/gene_essentiality_etfl_fast_1366.csv', index_col = 0)
# etfl_pred = pd.read_csv('outputs/gene_essentiality_vETFL_infer.csv', index_col = 0)
#etfl_pred = pd.read_csv('outputs/gene_essentiality_vETFL.csv', index_col = 0)
etfl_pred = pd.read_csv('outputs/gene_essentiality_vETFL_gpr.csv', index_col = 0)

classify = lambda df,x: df>x

df1 = pd.DataFrame(index = non_ess_genes['genes'])
df1['has_growth'] = True
df2 = pd.DataFrame(index = essential_genes['genes'])
df2['has_growth'] = False

ground_truth = pd.concat([df1, df2], axis=0, ignore_index=False)

common_genes = set(ground_truth.index).intersection(etfl_pred.index)

# mcc = dict()

# for x in linspace(0,etfl_pred.max().iloc[0],20):
#
#     y_true = ground_truth.loc[common_genes]
#     y_pred = classify(etfl_pred.loc[common_genes], x)
#
#     mcc[x] = matthews_corrcoef(y_true = y_true,
#                                y_pred = y_pred)

# threshold = etfl_pred.max().iloc[0] / 2

y_true = ground_truth.loc[common_genes]['has_growth']
y_pred = ~etfl_pred.loc[common_genes].isna()['0']
# y_pred = classify(etfl_pred.loc[common_genes], threshold)

mcc = matthews_corrcoef(y_true = y_true,
                           y_pred = y_pred)

cm = confusion_matrix(y_true = y_true,
                      y_pred = y_pred,
                      labels = [False, True])

# CM is Essential | Non essential <=> Has no growth || Has growth

calc_mcc = lambda tp,tn,fp,fn: (tp*tn - fp*fn)/ \
                               np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))


# From sklearn doc on confusion matrix
# Thus in binary classification, the count of true negatives is C[0,0],
# false negatives is C[1,0], true positives is C[1,1] and false positives is C[0,1].
# Predictions are columns, hence we need to transpose to have predictions as rows
cm = cm.T

print('iJO1366 Matthews correlation coefficient:',
      calc_mcc(len(ijo_tp),len(ijo_tn),len(ijo_fp),len(ijo_fn)))
print('iJO1366 Confusion matrix:')
print(np.array([[len(ijo_tn),len(ijo_fn)],[len(ijo_fp),len(ijo_tp)]]))
print('ETFL Matthews correlation coefficient:', mcc)
print('ETFL Confusion matrix:')
print(cm)

# Data for SI 6

etfl_tp  = y_true.loc[( y_true &  y_pred)]
etfl_tn  = y_true.loc[(~y_true & ~y_pred)]
etfl_fp  = y_true.loc[(~y_true &  y_pred)]
etfl_fn  = y_true.loc[( y_true & ~y_pred)]

etfl_classif = y_true.copy()
etfl_classif[etfl_tp.index] = 'tp'
etfl_classif[etfl_tn.index] = 'tn'
etfl_classif[etfl_fp.index] = 'fp'
etfl_classif[etfl_fn.index] = 'fn'

ijo_classif = y_true.copy()
ijo_classif[ijo_tp['genes']] = 'tp'
ijo_classif[ijo_tn['genes']] = 'tn'
ijo_classif[ijo_fp['genes']] = 'fp'
ijo_classif[ijo_fn['genes']] = 'fn'

from etfl.io.json import load_json_model
from etfl.core.reactions import EnzymaticReaction
m = load_json_model('models/SlackModel '
                            'iJO1366_vETFL_431_enz_128_bins__20190701_082518.json',
                    solver='cplex')

def enzrxn_to_gpr(rxn):
    return ' OR '.join([' AND '.join(['{}*{}'.format(v,k)
                                    for v,k in isozyme.composition.items()])
                                for isozyme in rxn.enzymes])
def enzrxn_to_iso(rxn):
    return ' OR '.join([isozyme.id for isozyme in rxn.enzymes])

list_of_offenders = list()
for gene_id in etfl_classif.index:
    if etfl_classif[gene_id] != ijo_classif[gene_id]:
        rxns = [x for x in m.reactions if any(gene_id == g.id for g in x.genes)]
        for rxn in rxns:
            ijo_GPR = rxn.gene_reaction_rule
            if isinstance(rxn, EnzymaticReaction):
                etfl_GPR = enzrxn_to_gpr(rxn)
                etfl_isoz = enzrxn_to_iso(rxn)
            else:
                etfl_GPR = 'no enzyme'
                etfl_isoz = ''
            record = [gene_id, rxn.id, rxn.subsystem,
                      ijo_classif[gene_id], etfl_classif[gene_id],
                      ijo_GPR, etfl_GPR, etfl_isoz]
            list_of_offenders.append(record)

df = pd.DataFrame(list_of_offenders,
                  columns = ['gene','reaction','subsystem',
                             'iJO_essentiality','vETFL_essentiality',
                             'iJO_GPR','ETFL_equivalent_GPR','vETFL_isozymes'])

df.to_excel('outputs/SI6_gene_ko_diffs.xlsx')