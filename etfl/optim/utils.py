# -*- coding: utf-8 -*-
"""
.. module:: etfl
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

Optimisation utilities
"""
from pytfa.optim.variables import ReactionVariable, MetaboliteVariable
from .variables import EnzymeVariable, GeneVariable, ModelVariable, \
    GrowthActivation, BinaryActivator
from pytfa.optim.constraints import ReactionConstraint, MetaboliteConstraint
from pytfa.optim.utils import get_all_subclasses
from .constraints import EnzymeConstraint, GeneConstraint, ModelConstraint
from collections import namedtuple


try:
    from gurobipy import GRB
except ModuleNotFoundError:
    pass


def make_subclasses_dict(cls):
    """
    Return a dictionary of the subclasses inheriting from the argument class.
    Keys are String names of the classes, values the actual classes.

    :param cls:
    :return:
    """
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
    """
    If the solver of the model whose integers to fix has Gurobi as a solver,
    use the built-in method

    :param model: A model with a Gurobi backend
    :return:
    """
    new = model.copy()
    fixed = model.solver.problem.fixed()
    new.solver.problem = fixed
    return new

def _generic_fix_integers(model):
    """
    Fix the integers of a model to its solution, and removes the variables.

    :param model:
    :return:
    """
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

        # TODO make fatser, using cons.change_expr and ad-hoc subs dicts
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
    """
    Rebuilds a variable from a classname and link it to the model

    :param classname:
    :param model:
    :param this_id:
    :param lb:
    :param ub:
    :param queue:
    :return:
    """
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
    """
    Rebuilds a constraint from a classname and link it to the model

    :param classname:
    :param model:
    :param this_id:
    :param new_expr:
    :param lb:
    :param ub:
    :param queue:
    :return:
    """
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
    """
    Check if the model uses Gurobi as a solver

    :param model:
    :return:
    """
    return model.problem.__name__ == 'optlang.gurobi_interface'


def fix_growth(model, solution = None):
    """
    Set the growth integers to their fixed values from a solution. If no
    solution is provided, the model's latest solution is used.
    The growth can be released using the function
    :func:`etfl.optim.utils.release_growth`

    :param model:
    :param solution:
    :return:
    """

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
    """
    Helper function. if solution is None, attempts to get it from the model.

    :param model:
    :param solution:
    :return:
    """
    if solution is None:
        try:
            solution = model.solution
        except AttributeError:
            raise AttributeError('If not providing a solution object, please '
                                 'provide a model with an embedded solution '
                                 '(call model.solve())')
    return solution


def release_growth(model):
    """
    After growth has been fixed by :func:`etfl.optim.utils.fix_growth`,
    it can be released using this function.

    :param model:
    :return:
    """

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
    """
    Gives a warm start to the model.
    Release it with :func:`etfl.optim.utils.release_warm_start`.

    :param model:
    :param solution:
    :return:
    """
    solution = check_solution(model, solution)

    if is_gurobi(model):
        for the_var in model.variables:
            if the_var.type == 'binary':
                the_var._internal_variable.Start = solution.raw[the_var.name]
    else:
        raise NotImplementedError('Solver not supported: ' + model.problem.__name__)


def release_warm_start(model):
    """
    Releases the warm start provided by
    :func:`etfl.optim.utils.apply_warm_start`.

    :param model:
    :return:
    """

    if is_gurobi(model):
        for the_var in model.variables:
            if the_var.type == 'binary':
                the_var._internal_variable.Start = GRB.UNDEFINED
    else:
        raise NotImplementedError('Solver not supported: ' + model.problem.__name__)


def get_active_growth_bounds(model):
    """
    Returns the growth bound closest to the growth flux calculated at the
    last solution.

    :param model:
    :return:
    """
    mu = model.growth_reaction.flux
    difflist = [abs(mu - x[0]) for x in model.mu_bins]
    min_diff = min(difflist)
    min_ix = difflist.index(min_diff)

    mu_i, (mu_lb, mu_ub) = model.mu_bins[min_ix]

    return mu_i, mu_lb, mu_ub


def safe_optim(model):
    """
    Catches *any* exception that can happen during solving, and logs it.
    Useful if you solve many problems in a sequence and some of then are
    infeasible.
    **Be careful** : This wil catch literally **any** Exception.

    :param model:
    :return:
    """
    try:
        out = model.optimize()
    except Exception as e:
        import numpy as np
        model.logger.warning('Exception occurred during solving: {}. - '
                             'Solver status: {}'.format(str(e), model.solver.status))
        out = DefaultSol
        out.objective_value = np.nan
    return  out