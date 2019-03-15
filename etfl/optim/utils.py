# -*- coding: utf-8 -*-
"""
.. module:: etfl
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: ETFL team

Optimisation utilities
"""
from pytfa.optim.variables import ReactionVariable, MetaboliteVariable
from .variables import EnzymeVariable, GeneVariable, ModelVariable, \
    GrowthActivation, BinaryActivator
from pytfa.optim.constraints import ReactionConstraint, MetaboliteConstraint
from .constraints import EnzymeConstraint, GeneConstraint, ModelConstraint
from collections import namedtuple

try:
    from gurobipy import GRB
except ModuleNotFoundError:
    pass


def get_all_subclasses(cls):
    all_subclasses = []

    for subclass in cls.__subclasses__():
        all_subclasses.append(subclass)
        all_subclasses.extend(get_all_subclasses(subclass))

    return all_subclasses

def make_subclasses_dict(cls):
    the_dict = {x.__name__:x for x in get_all_subclasses(cls)}
    the_dict[cls.__name__] = cls
    return the_dict

REACTION_VARIABLE_SUBCLASSES    = make_subclasses_dict(ReactionVariable)
REACTION_CONSTRAINT_SUBCLASSES  = make_subclasses_dict(ReactionConstraint)
METABOLITE_VARIABLE_SUBCLASSES  = make_subclasses_dict(MetaboliteVariable)
METABOLITE_CONSTRAINT_SUBCLASSES= make_subclasses_dict(MetaboliteConstraint)
ENZYME_VARIABLE_SUBCLASSES      = make_subclasses_dict(EnzymeVariable)
ENZYME_CONSTRAINT_SUBCLASSES    = make_subclasses_dict(EnzymeConstraint)
GENE_VARIABLE_SUBCLASSES        = make_subclasses_dict(GeneVariable)
GENE_CONSTRAINT_SUBCLASSES      = make_subclasses_dict(GeneConstraint)
MODEL_VARIABLE_SUBCLASSES       = make_subclasses_dict(ModelVariable)
MODEL_CONSTRAINT_SUBCLASSES     = make_subclasses_dict(ModelConstraint)

INTEGER_VARIABLE_TYPES = ('binary','integer')


def fix_integers(model):
    """
    Fixes all integer and binary variables of a model, to make it sample-able
    :param model:
    :return:
    """

    if model.problem.__name__ == 'optlang.gurobi_interface':
        model.logger.info('Gurobi-based model detected - using  Gurobi\'s .'
                          'fixed() method')
        return _gurobi_fix_integers(model)
    else:
        return _generic_fix_integers(model)

def _gurobi_fix_integers(model):
    new = model.copy()
    fixed = model.solver.problem.fixed()
    new.solver.problem = fixed
    return new

def _generic_fix_integers(model):
    continuous_model = model.copy()
    continuous_model.name = model.name + ' - continuous'

    integer_variables = set()

    constraints_with_integer_variables = []

    if not hasattr(model, 'solution'):
        model.logger.info('Model has no solution to fix the integers, calculating one')
        model.optimize()

    # We go through all the constraint descriptors and check if at least one of
    # their variables is in the integer variable list
    for this_cons in continuous_model._cons_dict.values():
        has_integer_variable = False
        for this_var in this_cons.constraint.variables:
            if this_var.type in INTEGER_VARIABLE_TYPES:
                has_integer_variable += True
                this_var_descriptor = this_var.name
                integer_variables.add(this_var_descriptor)
        constraints_with_integer_variables.append(this_cons.name)

    int_dict = {continuous_model.variables[x]: model.solution.x_dict[x]
                for x in integer_variables}

    for this_cons_name in constraints_with_integer_variables:
        this_cons = model._cons_dict[this_cons_name]
        new_expr = this_cons.expr.subs(int_dict)
        kind = type(this_cons)
        ub = this_cons.constraint.ub
        lb = this_cons.constraint.lb
        the_id = this_cons.id

        continuous_model.remove_constraint(this_cons)
        rebuild_constraint(classname=kind.__name__,
                           model=continuous_model,
                           this_id=the_id,
                           new_expr=new_expr,
                           lb=lb,
                           ub=ub)

    for this_var in integer_variables:
        # This_var is an InterfaceVariable object, we want the GenericVariable
        # it belongs to
        the_generic_var = continuous_model._var_dict[this_var.name]
        continuous_model.remove_variable(the_generic_var)

    continuous_model._push_queue()
    continuous_model.solver.update()
    # This will update the values =
    print('Is the cobra_model still integer ? {}'     \
          .format(continuous_model.solver.is_integer))

    return continuous_model

def rebuild_variable(classname, model, this_id, lb, ub, queue=True):
    if classname in REACTION_VARIABLE_SUBCLASSES:
        hook = model.reactions.get_by_id(this_id)
        this_class = REACTION_VARIABLE_SUBCLASSES[classname]
        nv = model.add_variable(kind=this_class,
                                hook=hook,
                                ub=ub,
                                lb=lb,
                                queue=queue)

    elif classname in METABOLITE_VARIABLE_SUBCLASSES:
        hook = model.metabolites.get_by_id(this_id)
        this_class = METABOLITE_VARIABLE_SUBCLASSES[classname]
        nv = model.add_variable(kind=this_class,
                                hook=hook,
                                ub=ub,
                                lb=lb,
                                queue=queue)

    elif classname in ENZYME_VARIABLE_SUBCLASSES:
        hook = model.enzymes.get_by_id(this_id)
        this_class = ENZYME_VARIABLE_SUBCLASSES[classname]
        nv = model.add_variable(kind=this_class,
                                hook=hook,
                                ub=ub,
                                lb=lb,
                                queue=queue)

    elif classname in GENE_VARIABLE_SUBCLASSES:
        hook = model.genes.get_by_id(this_id)
        this_class = GENE_VARIABLE_SUBCLASSES[classname]
        nv = model.add_variable(kind=this_class,
                                hook=hook,
                                ub=ub,
                                lb=lb,
                                queue=queue)

    elif classname in MODEL_VARIABLE_SUBCLASSES:
        hook = model
        this_class = MODEL_VARIABLE_SUBCLASSES[classname]
        nv = model.add_variable(kind=this_class,
                                hook=hook,
                                id_=this_id,
                                ub=ub,
                                lb=lb,
                                queue=queue)

    else:
        raise TypeError(
            'Class {} serialization not handled yet' \
                .format(classname))

    return nv


def rebuild_constraint(classname, model, this_id, new_expr, lb, ub, queue=True):
    if classname in REACTION_CONSTRAINT_SUBCLASSES:
        hook = model.reactions.get_by_id(this_id)
        this_class = REACTION_CONSTRAINT_SUBCLASSES[classname]
        nc = model.add_constraint(kind=this_class, hook=hook,
                                  expr=new_expr,
                                  ub=ub,
                                  lb=lb,
                                  queue=queue)

    elif classname in METABOLITE_CONSTRAINT_SUBCLASSES:
        hook = model.metabolites.get_by_id(this_id)
        this_class = METABOLITE_CONSTRAINT_SUBCLASSES[classname]
        nc = model.add_constraint(kind=this_class, hook=hook,
                                  expr=new_expr,
                                  ub=ub,
                                  lb=lb,
                                  queue=queue)

    elif classname in ENZYME_CONSTRAINT_SUBCLASSES:
        hook = model.enzymes.get_by_id(this_id)
        this_class = ENZYME_CONSTRAINT_SUBCLASSES[classname]
        nc = model.add_constraint(kind=this_class, hook=hook,
                                  expr=new_expr,
                                  ub=ub,
                                  lb=lb,
                                  queue=queue)

    elif classname in GENE_CONSTRAINT_SUBCLASSES:
        hook = model.genes.get_by_id(this_id)
        this_class = GENE_CONSTRAINT_SUBCLASSES[classname]
        nc = model.add_constraint(kind=this_class, hook=hook,
                                  expr=new_expr,
                                  ub=ub,
                                  lb=lb,
                                  queue=queue)

    elif classname in MODEL_CONSTRAINT_SUBCLASSES:
        hook = model
        this_class = MODEL_CONSTRAINT_SUBCLASSES[classname]
        nc = model.add_constraint(kind=this_class, hook=hook,
                                  expr=new_expr, id_=this_id,
                                  ub=ub,
                                  lb=lb,
                                  queue=queue)
    else:
        raise TypeError('Class {} serialization not handled yet' \
                        .format(classname))

    return nc


DefaultSol = namedtuple('DefaultSol', field_names='objective_value')


def is_gurobi(model):
    return model.problem.__name__ == 'optlang.gurobi_interface'


def fix_growth(model, solution = None):

    solution = check_solution(model, solution)

    mu_variables = model.get_variables_of_type(GrowthActivation)
    interp_variables = model.get_variables_of_type(BinaryActivator)

    vars_to_fix = list(mu_variables) + list(interp_variables)

    gurobi_hints = is_gurobi(model)
    if gurobi_hints:
        model.logger.info('Gurobi-based model detected - using  Gurobi hints')

    for the_var in vars_to_fix:
        value = solution.raw[the_var.name]
        try:
            the_var.variable.lb = int(value)
            the_var.variable.ub = int(value)
        except ValueError:
            # Happens if lb>ub during assignment
            the_var.variable.ub = int(value)
            the_var.variable.lb = int(value)

        if gurobi_hints:
            the_var.variable._internal_variable.VarHintVal = value
            the_var.variable._internal_variable.VarHintPri = 5


def check_solution(model, solution):
    if solution is None:
        try:
            solution = model.solution
        except AttributeError:
            raise AttributeError('If not providing a solution object, please '
                                 'provide a model with an embedded solution '
                                 '(call model.solve())')
    return solution


def release_growth(model):

    mu_variables = model.get_variables_of_type(GrowthActivation)
    interp_variables = model.get_variables_of_type(BinaryActivator)

    vars_to_fix = list(mu_variables) + list(interp_variables)

    gurobi_hints = is_gurobi(model)
    for the_var in vars_to_fix:
        the_var.variable.lb = 0
        the_var.variable.ub = 1

        if gurobi_hints:
            the_var.variable._internal_variable.VarHintVal = GRB.UNDEFINED
            the_var.variable._internal_variable.VarHintPri = 0


def apply_warm_start(model, solution):
    solution = check_solution(model, solution)

    for the_var in model.variables:
        if the_var.type == 'binary':
            the_var._internal_variable.Start = solution.raw[the_var.name]


def release_warm_start(model):

    for the_var in model.variables:
        if the_var.type == 'binary':
            the_var._internal_variable.Start = GRB.UNDEFINED


def get_active_growth_bounds(model):
    mu = model.growth_reaction.flux
    difflist = [abs(mu - x[0]) for x in model.mu_bins]
    min_diff = min(difflist)
    min_ix = difflist.index(min_diff)

    mu_i, (mu_lb, mu_ub) = model.mu_bins[min_ix]

    return mu_i, mu_lb, mu_ub


def safe_optim(model):
    try:
        out = model.optimize()
    except Exception:
        model.logger.warning('Solver status: {}'.format(model.solver.status))
        out = DefaultSol
        out.objective_value = np.nan
    return  out