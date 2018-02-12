import cobra
import numpy as np
from copy import deepcopy

from therme.therme.optim.constraints import CatalyticConstraint
from therme.therme.optim.variables import CatalyticActivator
from ..optim.constraints import BackwardDirectionCoupling
from ..optim.variables import BackwardUseVariable


def find_translation_gaps(model):
    tr_gaps = dict()
    for tr in model.translation_constraint:
        this_dict = {x: x.primal for x in tr.expr.free_symbols}
        this_dict['gap'] = tr.expr.evalf(subs=this_dict)
        tr_gaps[tr.name] = this_dict

    return tr_gaps


def find_essentials_from(model, met_dict):
    solvals = dict()
    for x, r in met_dict.items():
        previous_lb = r.lower_bound
        r.lower_bound = 0
        sol = model.optimize()
        solvals[x] = sol.f
        r.lower_bound = previous_lb
    return solvals


def get_model_argument(args, kwargs, arg_index = 0):
    try:
        model = kwargs['model']
    except KeyError:
        model = args[arg_index]

    return model

def save_objective_function(fun):
    def wrapper(*args, **kwargs):

        model = get_model_argument(args, kwargs)

        initial_obj = model.objective.expression

        out = fun(*args, **kwargs)

        model.objective = initial_obj

        return out
    return wrapper

def save_growth_bounds(fun):
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

    relaxation = deepcopy(model)

    relaxation.objective = relaxation.growth_reaction

    new_objective = 0

    for the_constr in relaxation.get_constraints_of_type(CatalyticConstraint).copy():
        activator = relaxation.add_variable(kind = CatalyticActivator,
                                               hook = the_constr.reaction,
                                               )
        new_expr = the_constr.constraint.expression - (1-activator) * relaxation.big_M

        lb = the_constr.constraint.lb
        ub = the_constr.constraint.ub

        relaxation.remove_constraint(the_constr)
        relaxation.add_constraint(kind = CatalyticConstraint,
                                     hook = the_constr.reaction,
                                     expr = new_expr,
                                     lb = lb,
                                     ub = ub)

        new_objective += activator

    relaxation.repair()

    relaxation.growth_reaction.lower_bound = min_growth

    relaxation.objective = new_objective

    relaxation.optimize()

    activators = relaxation.get_variables_of_type(CatalyticActivator)
    activator_states = relaxation.solution.x_dict[activators.list_attr('name')]

    for the_var in activators.list_attr("variable"):
        the_var.lb = int(relaxation.solution.x_dict[the_var.name])
        the_var.ub = int(relaxation.solution.x_dict[the_var.name])

    relaxed_model = deepcopy(model)

    for act in activators:
        if np.isclose(activator_states[act.name],0):
            the_cons = relaxed_model.catalytic_constraint.get_by_id(act.id)
            relaxed_model.remove_constraint(the_cons)

    return activator_states, relaxed_model, relaxation

