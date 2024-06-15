from copy import deepcopy

from pytfa.optim.utils import symbol_sum

import numpy as np

from optlang.symbolics import Zero, add

import cobra.util.solver as sutil
from cobra.util.solver import set_objective
from pytfa.optim.variables import ReactionVariable
from pytfa.optim.constraints import ReactionConstraint

class DeviationVariable(ReactionVariable):
    '''Set of variables defined to keep the difference between wild-type and knocke-out
    fluxes.'''
    
    prefix = 'DV_'
    
class HelperVariable(ReactionVariable):
    '''Set of variables defined to change absolute value to linear terms in 
    objective function.'''
    
    prefix = 'HV_'
    
class DeviationConstraint(ReactionConstraint):
    '''Set of constraints to determine the value of DeviationVariables'''
    
    prefix = 'DVC_'
    
class HelperConstraint1(ReactionConstraint):
    '''Set of constraints defined to change absolute value to linear terms in 
    objective function.'''
    
    prefix = 'HC1_'
    
class HelperConstraint2(ReactionConstraint):
    '''Set of constraints defined to change absolute value to linear terms in 
    objective function.'''
    
    prefix = 'HC2_'


def moma(model_perturbed, solution, linear=False, keep_model=True, scaled=False):
    '''Compute a single solution based on (linear) MOMA.

    Compute a new flux distribution that is at a minimal distance to a
    previous reference solution. Minimization of metabolic adjustment (MOMA) is
    generally used to assess the impact
    of knock-outs. Thus the typical usage is to provide a wildtype flux
    distribution as reference and a model_bkup in knock-out state.
    
    Parameters
    ----------
    model_perturbed : cobra.model_bkup, Thermo.model_bkup
        The model_bkup to add MOMA constraints and objective to.
    solution : cobra.Solution, optional
        A previous solution to use as a reference. If no solution is given,
        one will be computed using pFBA.
    linear : bool, optional
        Whether to use the linear MOMA formulation or not (default True).
    keep_model : bool, optional
        Whether to keep the current model_bkup or change it for MOMA formulation
    scaled : bool, optional
        Whether to scale the deviation for each flux. Scaling means to divide
        deviation by flux reference value in solution input.
        
    Notes
    -----
    In the original MOMA [1]_ specification one looks for the flux distribution
    of the deletion (v^d) closest to the fluxes without the deletion (v).
    In math this means:

    minimize \sum_i (v^d_i - v_i)^2
    s.t. Sv^d = 0
         lb_i <= v^d_i <= ub_i

    Here, we use a variable transformation v^t := v^d_i - v_i. Substituting
    and using the fact that Sv = 0 gives:

    minimize \sum_i (v^t_i)^2
    s.t. Sv^d = 0
         v^t = v^d_i - v_i
         lb_i <= v^d_i <= ub_i

    So basically we just re-center the flux space at the old solution and then
    find the flux distribution closest to the new zero (center). This is the
    same strategy as used in cameo.

    In the case of linear MOMA [2]_, we instead minimize \sum_i abs(v^t_i). The
    linear MOMA is typically significantly faster. Also quadratic MOMA tends
    to give flux distributions in which all fluxes deviate from the reference
    fluxes a little bit whereas linear MOMA tends to give flux distributions
    where the majority of fluxes are the same reference with few fluxes
    deviating a lot (typical effect of L2 norm vs L1 norm).

    '''
     
    if keep_model:
        model_bkup = deepcopy(model_perturbed)
    else:
        model_bkup = model_perturbed
        
    
    if not linear:
        model_bkup.solver = sutil.choose_solver(model_bkup, qp=True)
        
    # keeping the wild type growth rate
    setattr(model_bkup, 'original_growth', solution.objective_value)
    
    # resetting the objective function
    model_bkup.objective = model_bkup.problem.Objective(Zero, direction="min", sloppy=True)
    
    obj_vars = []
    for rxn in model_bkup.reactions:
        try:
            this_flux = solution.fluxes[rxn.id]
        except KeyError: # a reaction migh be added to the model_perturbed which is not in the wild-type
            this_flux = 0
        
        this_var = model_bkup.add_variable(kind = DeviationVariable, 
                                      hook = rxn, 
                                      lb = -1000, 
                                      ub = 1000)
        if scaled:
            if np.isclose(this_flux,0):
                sc_fact = 1/1e-5    # Basal flux
            else:
                sc_fact = 1/this_flux
            model_bkup.add_constraint(kind = DeviationConstraint, 
                                 hook = rxn, 
                                 expr = this_var \
                                 + model_bkup.variables[rxn.id]*sc_fact \
                                 - model_bkup.variables[rxn.reverse_id]*sc_fact,
                                 lb = 1, 
                                 ub = 1)
        else:
            model_bkup.add_constraint(kind = DeviationConstraint, 
                                 hook = rxn, 
                                 expr = this_var \
                                 + model_bkup.variables[rxn.id] \
                                 - model_bkup.variables[rxn.reverse_id],
                                 lb = this_flux, 
                                 ub = this_flux)

        if linear:
            # In this case we need to add one set of variables and two sets of 
            # constraints.
            helper_var = model_bkup.add_variable(kind = HelperVariable, 
                                                 hook = rxn, 
                                                 lb = 0, 
                                                 ub = 1000)
            model_bkup.add_constraint(kind = HelperConstraint1, 
                                      hook = rxn, 
                                      expr = this_var - helper_var,
                                      lb = -1000, 
                                      ub = 0)
            model_bkup.add_constraint(kind = HelperConstraint2, 
                                      hook = rxn, 
                                      expr = this_var + helper_var,
                                      lb = 0, 
                                      ub = 1000)
            obj_vars += [helper_var]
        else:
            obj_vars += [this_var]
    if linear:
        objective = symbol_sum([x for x in obj_vars])
    else:
        objective = add([x*x for x in obj_vars])
    
    set_objective(model_bkup, objective)
    moma_sol = model_bkup.optimize()
    
        
    return moma_sol