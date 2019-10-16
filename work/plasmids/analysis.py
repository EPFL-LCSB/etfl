from etfl.optim.utils import check_solution
from etfl.optim.variables import mRNAVariable, EnzymeVariable
from etfl.analysis.utils import enzymes_to_peptides_conc
from etfl.analysis.summary import print_standard_sol

from pytfa.io.viz import export_reactions_for_escher

import pandas as pd

from collections import OrderedDict


def chebyshev(model, chebyshev_include, inplace=True, big_m=1000):
    from etfl.optim.variables import mRNAVariable,EnzymeVariable
    from etfl.optim.constraints import ForwardCatalyticConstraint,BackwardCatalyticConstraint
    from etfl.analysis.dynamic import chebyshev_center, compute_center

    chebyshev_variables = model.get_variables_of_type(mRNAVariable)
    chebyshev_variables += model.get_variables_of_type(EnzymeVariable)

    initial_solution = model.optimize()

    the_obj = model.objective

    chebyshev_center(model, chebyshev_variables,
                 inplace=inplace,
                 big_m=big_m,
                 include=chebyshev_include)

    # Revert to the previous objective, as chebyshev_center sets it to maximize
    # the radius
    model.objective = the_obj

    # We want to compute the center under the constraint of the growth at the
    # initial solution.
    chebyshev_sol = compute_center(model, provided_solution=initial_solution,
                                   revert_changes=False)
    print_standard_sol(model)
    return chebyshev_sol


def export_mrna(model, tag, solution=None):
    solution = check_solution(model, solution)

    mrna_var_ids = model.get_variables_of_type(mRNAVariable).list_attr('name')

    data = solution.raw
    mrna_data = data[mrna_var_ids]
    mrna_data.name = tag

    filename = 'outputs/' + tag + '_mrna' + '.csv'

    mrna_data.to_csv(filename, header=True)
    return mrna_data

def export_peptides(model, tag, solution=None):
    solution = check_solution(model, solution)

    enz_var_ids = model.get_variables_of_type(EnzymeVariable).list_attr('name')

    data = solution.raw
    enz_data = data[enz_var_ids]
    enz_data.name = tag

    pep_data = enzymes_to_peptides_conc(model=model, enzyme_conc=enz_data)
    pep_data = pd.DataFrame.from_dict(pep_data, orient='index')

    filename = 'outputs/' + tag + '_pep' + '.csv'

    pep_data.to_csv(filename, header=True)
    return pep_data


def export_fluxes(model, tag, solution=None):
    solution = check_solution(model, solution)

    filename = 'outputs/' + tag + '_rxn' + '.csv'

    export_reactions_for_escher(model, solution.raw, filename)
    return solution.fluxes


def binding_constraints(model, tag, epsilon, solution=None):
    from etfl.optim.utils import get_binding_constraints
    # solution = check_solution(model, solution)
    model.optimize()

    filename = 'outputs/' + tag + '_binding_cons' + '.csv'

    binding_constraints = get_binding_constraints(model, epsilon)

    name2bind = OrderedDict()

    for cons_kind, cons_names in binding_constraints.items():
        for the_name in cons_names:
            name2bind[the_name] = cons_kind

    df = pd.DataFrame.from_dict(name2bind, orient='index')
    df.columns = ['name','type']

    df.to_csv(filename, header=True)
    return df
