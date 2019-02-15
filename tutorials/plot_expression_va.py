import os
import pandas as pd
import re
import numpy as np
import bokeh.plotting as bp
from bokeh.models.glyphs import HBar, VBar
from bokeh.models import ColumnDataSource, DataRange1d, LinearAxis, \
    Grid,CategoricalTicker, FuncTickFormatter, HoverTool
from bokeh.io import export_svgs

outputs_folder = './outputs/'
plots_folder = './plots/'

verbose_kinds = {
    'mrna':'mRNA',
    'rxns':'Reaction',
    'enz': 'Enzyme',
    'pep': 'Peptide',
}

VA_REGEX = re.compile(r'iJO_(T[01]E[01]N[01])_low?_hi_\-?\d+\.?\d*_({})'
                        .format('|'.join(verbose_kinds.keys())))

# DEFAULT_COLOR = "#b3de69" # Light Green
DEFAULT_COLOR = "#2579b6" # Mid blue

def get_va_files():

    files = dict()

    for filename in os.listdir(outputs_folder):
        re_search = VA_REGEX.search(filename)
        if re_search and filename.endswith('.csv'):
            tag, kind = re_search.groups(1)
            files[filename] = (tag, kind)

    return files

def plot_va(filename, tag, kind, color = DEFAULT_COLOR, orient = 'horizontal'):
    bp.curdoc().clear()

    title = verbose_kinds[kind] + ' variability analysis for {} iJO1366 model'\
                                    .format(tag)
    data = pd.read_csv(os.path.join(outputs_folder, filename), index_col = 0)

    if not data.columns[0] in ['minimum','maximum']:
        data.columns = ['minimum','maximum']

    data += 1e-9 # Resolution of the solver
    f = lambda x: np.sqrt(x[0]*x[1])
    # data['score'] = data.mean(axis=1)
    data['score'] = data.apply(f, axis=1)
    data.sort_values(by='score', ascending = False, inplace = True)
    data['y'] = range(len(data))
    data['name'] = data.index

    source = ColumnDataSource(data)

    xdr = DataRange1d()
    ydr = DataRange1d()

    _tools_to_show = 'box_zoom,pan,save,hover,reset,tap,wheel_zoom'


    if orient == 'vertical':
        p1 = bp.figure( title=title, x_range=xdr, y_range=ydr,
                        x_axis_type = 'log',
                        plot_width=600,
                        plot_height=1000,
                        tools=_tools_to_show,
                        # h_symmetry=False, v_symmetry=False,
                        min_border=0)
        glyph = HBar(y="y", right="maximum", left="minimum", height=0.9,
                     fill_color=color, fill_alpha=1,
                     line_color = None)
        p1.add_glyph(source, glyph)

        p1.circle(x='score', y='y', fill_color='white', line_color=color,
                  source=source)
        axis1 = p1.xaxis
        axis2 = p1.yaxis
    elif orient == 'horizontal':
        p1 = bp.figure(title=title, x_range=ydr, y_range=xdr,
                       y_axis_type='log',
                       plot_width=1000,
                       plot_height=600,
                       tools=_tools_to_show,
                       # h_symmetry=False, v_symmetry=False,
                       min_border=0)

        glyph = VBar(x="y", top="maximum", bottom="minimum", width=0.9,
                     fill_color=color, fill_alpha=1,
                     line_color=None)
        p1.add_glyph(source, glyph)

        p1.circle(y='score', x='y', fill_color='white', line_color=color,
                  source=source)
        axis1 = p1.yaxis
        axis2 = p1.xaxis
    else:
        raise ValueError("orient should be 'vertical' or 'horizontal'")

    # Fix ticks
    label_dict = {}
    for i, s in enumerate(data.index):
        label_dict[i] = s
    axis2.formatter = FuncTickFormatter(code="""
                                            var labels = %s;
                                            return labels[tick];
                                        """ % label_dict)
    axis1.axis_label = '[{}]'.format(verbose_kinds[kind])

    # p1.yaxis.ticker = [x for x in range(len(data))]

    hover = p1.select(dict(type=HoverTool))
    hover.tooltips = [(verbose_kinds[kind], "@name"),
                      ("min", "@minimum"),
                      ("max", "@maximum"),
                      ]
    hover.mode = 'mouse'


    path = os.path.join(plots_folder, 'va_' + filename)
    bp.output_file(path + '.html')
    bp.show(p1)
    p1.output_backend='svg'
    export_svgs(p1, filename=path+'.svg')
    bp.curdoc().clear()

    return data


if __name__ == '__main__':
    the_files = get_va_files()

    for filename, (tag, kind) in the_files.items():
        data = plot_va(filename, tag, kind)