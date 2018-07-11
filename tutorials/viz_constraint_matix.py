import pandas as pd
import bokeh.plotting as bp
from bokeh.models import HoverTool, ColumnDataSource, Legend


"""
Adapted from
https://stackoverflow.com/a/40231779
"""

tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

def get_expr_coos(expr, var_indices):
    for i in range(expr.size()):
        dvar = expr.getVar(i)
        yield dvar, expr.getCoeff(i), var_indices[dvar]


def get_matrix_coo_gurobi(model):
    m = model.solver.problem

    dvars = m.getVars()
    constrs = m.getConstrs()
    var_indices = {v: i for i, v in enumerate(dvars)}
    for row_idx, constr in enumerate(constrs):
        for this_var, coeff, col_idx in get_expr_coos(m.getRow(constr), var_indices):
            vartype = this_var.VType
            varname = this_var.VarName
            yield row_idx, col_idx, coeff, varname, vartype, constr.ConstrName


def get_constraint_matrix_gurobi(model):

    dvars = m.getVars()
    constrs = m.getConstrs()

    obj_coeffs = m.getAttr('Obj', dvars)

    var_index = {v: i for i, v in enumerate(dvars)}
    constr_index= {c: i for i, c in enumerate(constrs)}


    nzs = pd.DataFrame(get_matrix_coo_gurobi(m),
                       columns=['row_idx',
                                'col_idx',
                                'coeff',
                                'var_name',
                                'var_type',
                                'cons_name'])

    return nzs

def add_color_group(nzs):

    groups = pd.Series(index = nzs.index)
    groups.fillna('0ther', inplace = True)

    # groups[nzs['name'].str.contains('transcription')] = 'transcription'
    # groups[nzs['name'].str.contains('translation')] = 'translation'
    #
    # groups[nzs['name'].str.startswith('MR_')] = 'mRNA'
    # groups[nzs['name'].str.startswith('EZ_')] = 'Enzyme'
    # groups[nzs['name'].str.startswith('GA_')] = 'Binary Growth coefficient'
    # groups[nzs['name'].str.startswith('LZ_')] = 'Linearization variable'

    groups[nzs['cons_name'].str.startswith('TR_')] = 'Transcription/translation'
    groups[nzs['cons_name'].str.startswith('LC_')] = 'Linearization'
    groups[nzs['cons_name'].str.startswith('MB_')] = 'mRNA balance'
    groups[nzs['cons_name'].str.startswith('MD_')] = 'mRNA degradation'
    groups[nzs['cons_name'].str.startswith('EB_')] = 'Enzyme balance'
    groups[nzs['cons_name'].str.startswith('ED_')] = 'Enzyme degradation'
    groups[nzs['cons_name'].str.startswith('FC_')] = 'Enzyme kcat'
    groups[nzs['cons_name'].str.startswith('BC_')] = 'Enzyme kcat'
    groups[nzs['cons_name'].str.startswith('TC_')] = 'Capacity'
    groups[nzs['cons_name'].str.startswith('EX_')] = 'Transcription/translation coupling'

    groups[nzs['cons_name'].str.startswith('b')] = 'Peptide'

    nzs['groups'] = groups

    return nzs


def plot_spy_matrix(model, destination):

    nzs = pd.DataFrame(get_matrix_coo_gurobi(model),
                       columns=['row_idx',
                                'col_idx',
                                'coeff',
                                'var_name',
                                'var_type',
                                'cons_name'])

    nzs = add_color_group(nzs)
    nzs['row_idx'] *= -1

    bp.curdoc().clear()

    hover = HoverTool(tooltips=[("constraint name", "@cons_name"),
                                ("variable name", "@var_name"),
                                ("coefficient", "@coeff"),
                                ("(x,y)", "(@col_idx, @row_idx)"),
                                ("type", "@var_type"), ])

    p = bp.figure(width = 1000)
    p.add_tools(hover)

    legend_it = []

    for e, (subgroup, this_data) in enumerate(nzs.groupby('groups')):

        source = ColumnDataSource.from_df(this_data)

        # if this_data['var_type'].iloc[0] in ['I', 'B']:
        #     plot_fun = p.square
        # else:
        #     plot_fun = p.circle

        plot_fun = p.square
        c = plot_fun('col_idx', 'row_idx',
                     color = tableau20[e],
                     source = source,
                     size = 2,
                     alpha=1)

        legend_it.append((subgroup,[c]))

    legend = Legend(items=legend_it, location=(0, 0))
    legend.click_policy = "mute"

    p.add_layout(legend, 'right')

    bp.output_file(destination)
    bp.show(p)

    return p