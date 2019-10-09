from etfl.optim.utils import check_solution
from etfl.optim.variables import mRNAVariable, EnzymeVariable


def chebyshev(model, chebyshev_include, inplace=True, big_m=1000):
    from etfl.optim.variables import mRNAVariable,EnzymeVariable
    from etfl.optim.constraints import ForwardCatalyticConstraint,BackwardCatalyticConstraint
    from etfl.analysis.dynamic import chebyshev_center, compute_center

    chebyshev_variables = model.get_variables_of_type(mRNAVariable)
    chebyshev_variables += model.get_variables_of_type(EnzymeVariable)

    chebyshev_center(model, chebyshev_variables,
                 inplace=inplace,
                 big_m=big_m,
                 include=chebyshev_include)

    # Revert to the previous objective, as chebyshev_center sets it to maximize
    # the radius
    model.objective = the_obj

    # We want to compute the center under the constraint of the growth at the
    # initial solution.
    chebyshev_sol = compute_center(model, provided_solution=initial_solution)
    print_standard_sol(model)


def export_mrna(model, tag, solution=None):
    solution = check_solution(solution)

    mrna_var_ids = model.get_variables_of_type(mRNAVariable).list_attr('id')

    data = solution.raw
    mrna_data = data[mrna_var_ids]
    mrna_data.name = tag

    filename = 'outputs/' + tag + 'mrna' + '.csv'

    mrna_data.to_csv(filename, header=True)