from etfl.optim.utils import check_solution
from etfl.optim.variables import mRNAVariable, EnzymeVariable, DNAVariable, RibosomeUsage, RNAPUsage
from etfl.analysis.utils import enzymes_to_peptides_conc
from etfl.analysis.summary import print_standard_sol

from pytfa.io.viz import export_reactions_for_escher

import pandas as pd

from collections import OrderedDict

def mkdir(directory):
    import os
    full_dir = 'outputs/' + directory
    if not os.path.exists(full_dir):
        os.makedirs(full_dir)

def chebyshev(model, chebyshev_include, inplace=True, big_m=1000):
    from etfl.optim.variables import mRNAVariable,EnzymeVariable, RibosomeUsage, RNAPUsage
    from etfl.optim.constraints import ForwardCatalyticConstraint,BackwardCatalyticConstraint
    from pytfa.analysis.chebyshev import chebyshev_transform, \
        get_cons_var_classes, get_variables
    from etfl.analysis.dynamic import compute_center

    if not inplace:
        new = model.copy()
    else:
        new = model

    new.objective = new.growth_reaction
    initial_solution = new.optimize()
    print('Initial solution:')
    print_standard_sol(new,solution=initial_solution)

    # List of radii to add
    radii = list()

    the_obj = new.objective

    for this_cons_type in chebyshev_include:
        if this_cons_type == 'CatalyticConstraint':
            # chebyshev_variables = new.get_variables_of_type(EnzymeVariable)
            # vars = get_variables(new, chebyshev_variables)
            include_list = get_cons_var_classes(new, ['ForwardCatalyticConstraint', 'BackwardCatalyticConstraint'],
                                                type = 'cons')
            this_r_id = 'catalytic'
            scaling_factor=1#172*3600
        elif this_cons_type =='ExpressionCoupling':
            # chebyshev_variables = new.get_variables_of_type(mRNAVariable)
            # chebyshev_variables = new.get_variables_of_type(RibosomeUsage)
            # vars = get_variables(new, chebyshev_variables)
            include_list = get_cons_var_classes(new, [this_cons_type], type='cons')
            this_r_id = 'translation'
            scaling_factor = 1#50
        elif this_cons_type =='RNAPAllocation':
            # chebyshev_variables = new.get_variables_of_type(DNAVariable)
            # chebyshev_variables = new.get_variables_of_type(RNAPUsage)
            # vars = get_variables(new, chebyshev_variables)
            include_list = get_cons_var_classes(new, [this_cons_type], type='cons')
            this_r_id = 'transcription'
            scaling_factor = 1
        elif this_cons_type =='SynthesisConstraint':
            # chebyshev_variables = new.get_variables_of_type(DNAVariable)
            # chebyshev_variables = new.get_variables_of_type(RNAPUsage)
            # vars = get_variables(new, chebyshev_variables)
            include_list = get_cons_var_classes(new, [this_cons_type], type='cons')
            this_r_id = 'synthesis'
            scaling_factor = 1
        else:
            include_list = get_cons_var_classes(new, [this_cons_type], type='cons')
            this_r_id = this_cons_type
            scaling_factor = 1


        r = chebyshev_transform(model=new,
                                vars=list(),
                                # vars=vars,
                                include_list=include_list,
                                exclude_list=list(),
                                radius_id=this_r_id, # if commented, all the radii will have
                                # the same name, hence it will be the same var
                                scaling_factor=scaling_factor)
        radii.append(r)

    regularisation =  0.001*new.enzymes.dummy_enzyme.scaled_concentration
    cheby_obj = sum(radii) #+ regularisation

    # We want to compute the center under the constraint of the growth at the
    # initial solution.
    chebyshev_sol = compute_center(new, provided_solution=initial_solution,
                                   objective=cheby_obj,
                                   revert_changes=False)
    print_standard_sol(new, solution=chebyshev_sol)
    return chebyshev_sol


def export_mrna(model, tag, solution=None):
    solution = check_solution(model, solution)

    mrna_var_ids = model.get_variables_of_type(mRNAVariable).list_attr('name')

    data = solution.raw
    mrna_data = data[mrna_var_ids]
    mrna_data.name = tag

    mkdir(tag)
    filename = 'outputs/' + tag + '/mrna' + '.csv'

    mrna_data.to_csv(filename, header=True)
    return mrna_data

def export_peptides(model, tag, solution=None):
    solution = check_solution(model, solution)

    enz_var_ids = model.get_variables_of_type(EnzymeVariable).list_attr('name')

    data = solution.raw
    enz_data = data[enz_var_ids]
    enz_data.columns = [tag]

    pep_data = enzymes_to_peptides_conc(model=model, enzyme_conc=enz_data)
    pep_data = pd.DataFrame.from_dict(pep_data, orient='index')

    mkdir(tag)
    filename = 'outputs/' + tag + '/pep' + '.csv'

    pep_data.to_csv(filename, header=True)
    return pep_data


def export_fluxes(model, tag, solution=None):
    solution = check_solution(model, solution)

    mkdir(tag)
    filename = 'outputs/' + tag + '/rxn' + '.csv'

    export_reactions_for_escher(model, solution.raw, filename)
    return solution.fluxes


def binding_constraints(model, tag, epsilon, solution=None):
    from etfl.optim.utils import get_binding_constraints
    # solution = check_solution(model, solution)
    model.optimize()

    mkdir(tag)
    filename = 'outputs/' + tag + '/binding_cons' + '.csv'

    binding_constraints = get_binding_constraints(model, epsilon)

    name2bind = OrderedDict()

    for cons_kind, cons_names in binding_constraints.items():
        for the_name in cons_names:
            name2bind[the_name] = cons_kind

    df = pd.DataFrame.from_dict(name2bind, orient='index')
    df.columns = ['type']
    df.index.name = 'name'

    df.to_csv(filename, header=True)
    return df


def export_slacks(model, tag, constraint_types, solution=None):
    from pytfa.analysis.chebyshev import get_cons_var_classes, ChebyshevRadius

    solution = check_solution(model, solution)
    include_list = get_cons_var_classes(model, constraint_types,
                                type = 'cons')
    cons_dict = {}

    if model.problem.__name__ == 'optlang.gurobi_interface':
        for cons_class in include_list:
            these_cons = model.get_constraints_of_type(cons_class)
            cons_dict.update({c.name:(c.constraint._internal_constraint.Slack,
                                      get_radius(c.expr),
                                      c.__class__.__name__)
                               for c in these_cons})

    mkdir(tag)
    filename = 'outputs/' + tag + '/slacks' + '.csv'

    df = pd.DataFrame.from_dict(cons_dict, orient='index')
    df.columns = ['slack','radius','type']
    df.index.name = 'name'

    df.to_csv(filename, header=True)
    return df

def get_radius(expr):
    """
    Returns the value of the sum of Chebyshev radii in a constraint expression

    :param cons:
    :return:
    """
    from pytfa.analysis.chebyshev import ChebyshevRadius
    variables = expr.free_symbols
    total = sum([v.primal for v in variables
                 if v.name.startswith(ChebyshevRadius.prefix)])
    return total