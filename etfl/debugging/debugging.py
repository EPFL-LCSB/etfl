import cobra
import numpy as np
from copy import deepcopy

from ..optim.constraints import CatalyticConstraint
from ..optim.variables import CatalyticActivator
from pytfa.optim.constraints import BackwardDirectionCoupling
from pytfa.optim.variables import BackwardUseVariable
from ..optim.constraints import ForwardCatalyticConstraint, BackwardCatalyticConstraint

def localize_exp(exp):
    """
    Takes an optlang expression, and replaces symbols (tied to variables) by
    their string names, to compare expressions of two different models

    :param exp:
    :type exp: :class:`optlang.symbolics.Expr`
    :return:
    """
    return exp.subs({x : x.name for x in exp.free_symbols})

def compare_expressions(exp1, exp2):
    """
    Check is the two given expressions are equal

    :param exp1:
    :type exp1: :class:`optlang.symbolics.Expr`
    :param exp2:
    :type exp2: :class:`optlang.symbolics.Expr`
    :return:
    """
    local_1 = localize_exp(exp1)
    local_2 = localize_exp(exp2)

    return local_1 == local_2

def find_different_constraints(model1, model2):
    """
    Given two models, find which expressions are different

    :param model1:
    :param model2:
    :return:
    """
    out = []
    for x in model1.constraints:
        y = model2.constraints.get(x.name)
        l1 = localize_exp(x.expression)
        l2 = localize_exp(x.expression)
        if l1 != l2:
            out += [x.name, l1, l2]
    return out

def find_translation_gaps(model):
    """
    For each translation constraint in the model, finds the value of each
    variable, and then evaluates the LHS of the constraint

    Constraints look like
    v_tsl - ktrans/L [RNAP_i] <= 0

    :param model:
    :return:
    """
    tr_gaps = dict()
    for tr in model.translation_constraint:
        this_dict = {x: x.primal for x in tr.expr.free_symbols}
        this_dict['gap'] = tr.expr.evalf(subs=this_dict)
        tr_gaps[tr.name] = this_dict

    return tr_gaps


def find_essentials_from(model, met_dict):
    """
    Given a dictionnary of {met_id:uptake_reaction}, checks the value of the
    objective function at optimality when the given uptake is closed.

    **Uptake reactions are expected to be aligned according to the consensus
    directionality for systems : met_e <=> []**

    :param model:
    :param met_dict:
    :return:
    """
    solvals = dict()
    for x, r in met_dict.items():
        previous_lb = r.lower_bound
        r.lower_bound = 0
        sol = model.optimize()
        solvals[x] = sol.f
        r.lower_bound = previous_lb
    return solvals


def get_model_argument(args, kwargs, arg_index = 0):
    """
    Utility function to get the model object from the arguments of a function

    :param args:
    :param kwargs:
    :param arg_index:
    :return:
    """
    try:
        model = kwargs['model']
    except KeyError:
        model = args[arg_index]

    return model

def save_objective_function(fun):
    """
    Decorator to restore the objective function after the execution of the
    decorated function.

    :param fun:
    :return:
    """
    def wrapper(*args, **kwargs):

        model = get_model_argument(args, kwargs)

        initial_obj = model.objective.expression

        out = fun(*args, **kwargs)

        model.objective = initial_obj

        return out
    return wrapper

def save_growth_bounds(fun):
    """
    Decorator to save the growth bound and restore them after the execution of
    the decorated function.

    :param fun:
    :return:
    """
    def wrapper(*args, **kwargs):

        model = get_model_argument(args, kwargs)

        initial_growth_lb = model.growth_reaction.lower_bound

        out = fun(*args, **kwargs)

        model.growth_reaction.lower_bound = initial_growth_lb

        return out

    return wrapper

@save_objective_function
@save_growth_bounds
def perform_iMM(model, uptake_dict, min_growth_coef=0.5, bigM=1000):
    """
    An implementation of the *in silico Minimal Media* methods, which uses MILP
    to find a minimum set of uptakes necessary to meet growth requirements

    See:

    Bioenergetics-based modeling of Plasmodium falciparum metabolism reveals
    its essential genes, nutritional requirements, and thermodynamic bottlenecks
    Chiappino-Pepe A, Tymoshenko S, Ataman M, Soldati-Favre D, Hatzimanikatis V
    (2017) PLOS Computational Biology 13(3): e1005397.
    https://doi.org/10.1371/journal.pcbi.1005397

    :param model:
    :type model: :class:`etfl.core.memodel.MEModel`
    :param uptake_dict: {met_id : <reaction object>}
    :param min_growth_coef: minimum fraction of optimal growth to be met
    :param bigM: a big-M value for the optimization problem
    :return:
    """
    use_vars = dict()
    use_cons = dict()


    for met, rxn in uptake_dict.items():
        BU_rxn = model.add_variable(kind=BackwardUseVariable, hook=rxn)

        # UR_rxn: R_rxn - M RU_rxn < 0
        R_rxn = rxn.reverse_variable
        CLHS = R_rxn - BU_rxn * bigM

        UC = model.add_constraint(kind=BackwardDirectionCoupling, hook=rxn,
                                  expr=CLHS, ub=0)

        use_vars[met] = BU_rxn
        use_cons[met] = UC

    initial_growth = model.growth_reaction.forward_variable.primal

    model.growth_reaction.lower_bound = min_growth_coef * initial_growth

    sum_use = sum(use_vars.values())

    model.objective = -1 * sum_use

    model.optimize()

    minimal_medium = {x: r.reverse_variable.primal for x, r in
                      uptake_dict.items()}

    for cons in use_cons.values():
        model.remove_constraint(cons)
    for var in use_vars.values():
        model.remove_variable(var)

    return minimal_medium


@save_objective_function
def check_production_of_mets(model, met_ids):
    """
    for each metabolite ID given, create a sink and maximize the production of
    the metabolite

    :param model:
    :type model: :class:`etfl.core.memodel.MEModel`
    :param met_ids:
    :return:
    """
    # Check that the FBA model can produce each metabolite in the given list
    met_production = dict()

    for met_id in met_ids:
        # # Free amino acids:
        met_c = model.metabolites.get_by_id(met_id)

        met_e = cobra.Metabolite(id=met_id + '_e', name=met_id, compartment='e')
        met_exchange = cobra.Reaction(id='EX_' + met_id,
                                      name=met_id + ' Exchange',
                                      lower_bound=-10, upper_bound=100)
        met_exchange.add_metabolites({met_e: 1, met_c: -1})

        met_drain = cobra.Reaction(id='DM_' + met_id, name=met_id + ' Drain',
                                   lower_bound=0, upper_bound=100)
        met_drain.add_metabolites({met_e: -1})

        model.add_reactions([met_exchange, met_drain])
        the_rxn = model.reactions.get_by_id(met_drain.id)
        model.objective = the_rxn

        met_production[met_id] = model.optimize().f

    return met_production


def relax_catalytic_constraints(model, min_growth):
    """
    Find a minimal set of catalytic constraints to relax to meet a minimum
    growth criterion

    :param model:
    :type model: :class:`etfl.core.memodel.MEModel`
    :param min_growth:
    :return:
    """

    relaxation = deepcopy(model)

    relaxation.objective = relaxation.growth_reaction

    new_objective = 0

    for the_constr in relaxation.forward_catalytic_constraint:
        activator = relaxation.add_variable(kind = CatalyticActivator,
                                               hook = the_constr.reaction,
                                               )
        
        new_expr = the_constr.constraint.expression - (1-activator) * relaxation.big_M

        lb = the_constr.constraint.lb
        ub = the_constr.constraint.ub

        relaxation.remove_constraint(the_constr)
        relaxation.add_constraint(kind = ForwardCatalyticConstraint,
                                     hook = the_constr.reaction,
                                     expr = new_expr,
                                     lb = lb,
                                     ub = ub)
        
        # we should consider the corresponding backward constraint
#        bkwd_id=the_constr.reaction.id
#        the_constr_bkwd = relaxation.backward_catalytic_constraint.get_by_id(bkwd_id)
#        new_expr_bkwd = the_constr_bkwd.constraint.expression - (1-activator) * relaxation.big_M
#        
#        lb_bkwd = the_constr_bkwd.constraint.lb
#        ub_bkwd = the_constr_bkwd.constraint.ub
#        
#        relaxation.remove_constraint(the_constr_bkwd)
#        relaxation.add_constraint(kind = BackwardCatalyticConstraint,
#                                     hook = the_constr_bkwd.reaction,
#                                     expr = new_expr_bkwd,
#                                     lb = lb_bkwd,
#                                     ub = ub_bkwd)

        new_objective += activator

    relaxation.repair()

    relaxation.growth_reaction.lower_bound = min_growth

    relaxation.objective = new_objective

    relaxation.optimize()

    activators = relaxation.get_variables_of_type(CatalyticActivator)
    activator_states=relaxation.solution.raw.loc[activators.list_attr('name')]

    for the_var in activators.list_attr("variable"):
        the_var.lb = int(relaxation.solution.raw.loc[the_var.name])
        the_var.ub = int(relaxation.solution.raw.loc[the_var.name])

    relaxed_model = deepcopy(model)

    for act in activators:
        if np.isclose(activator_states[act.name],0):
            the_cons = relaxed_model.forward_catalytic_constraint.get_by_id(act.id)
            relaxed_model.remove_constraint(the_cons)

    return activator_states, relaxed_model, relaxation

def relax_catalytic_constraints_bkwd(model, min_growth):
    """
    Find a minimal set of catalytic constraints to relax to meet a minimum
    growth criterion

    :param model:
    :type model: :class:`etfl.core.memodel.MEModel`
    :param min_growth:
    :return:
    """

    relaxation = deepcopy(model)

    relaxation.objective = relaxation.growth_reaction

    new_objective = 0

    for the_constr in relaxation.backward_catalytic_constraint:
        activator = relaxation.add_variable(kind = CatalyticActivator,
                                               hook = the_constr.reaction,
                                               )
        
        new_expr = the_constr.constraint.expression - (1-activator) * relaxation.big_M

        lb = the_constr.constraint.lb
        ub = the_constr.constraint.ub

        relaxation.remove_constraint(the_constr)
        relaxation.add_constraint(kind = BackwardCatalyticConstraint,
                                     hook = the_constr.reaction,
                                     expr = new_expr,
                                     lb = lb,
                                     ub = ub)
        
        # we should consider the corresponding backward constraint
#        bkwd_id=the_constr.reaction.id
#        the_constr_bkwd = relaxation.backward_catalytic_constraint.get_by_id(bkwd_id)
#        new_expr_bkwd = the_constr_bkwd.constraint.expression - (1-activator) * relaxation.big_M
#        
#        lb_bkwd = the_constr_bkwd.constraint.lb
#        ub_bkwd = the_constr_bkwd.constraint.ub
#        
#        relaxation.remove_constraint(the_constr_bkwd)
#        relaxation.add_constraint(kind = BackwardCatalyticConstraint,
#                                     hook = the_constr_bkwd.reaction,
#                                     expr = new_expr_bkwd,
#                                     lb = lb_bkwd,
#                                     ub = ub_bkwd)

        new_objective += activator

    relaxation.repair()

    relaxation.growth_reaction.lower_bound = min_growth

    relaxation.objective = new_objective

    relaxation.optimize()

    activators = relaxation.get_variables_of_type(CatalyticActivator)
    activator_states=relaxation.solution.raw.loc[activators.list_attr('name')]

    for the_var in activators.list_attr("variable"):
        the_var.lb = int(relaxation.solution.raw.loc[the_var.name])
        the_var.ub = int(relaxation.solution.raw.loc[the_var.name])

    relaxed_model = deepcopy(model)

    for act in activators:
        if np.isclose(activator_states[act.name],0):
            the_cons = relaxed_model.backward_catalytic_constraint.get_by_id(act.id)
            relaxed_model.remove_constraint(the_cons)

    return activator_states, relaxed_model, relaxation