from copy import deepcopy

from pytfa.optim.utils import symbol_sum

import numpy as np

from optlang.symbolics import Zero, add

import cobra.util.solver as sutil
from cobra.util.solver import set_objective
from etfl.optim.variables import EnzymeVariable
from etfl.optim.constraints import EnzymeConstraint

class DeviationVariable(EnzymeVariable):
    '''Set of variables defined to keep the difference between wild-type and knocke-out
    enzyme values.'''
    
    prefix = 'DV_'
    
class HelperVariable(EnzymeVariable):
    '''Set of variables defined to change absolute value to linear terms in 
    objective function.'''
    
    prefix = 'HV_'
    
class DeviationConstraint(EnzymeConstraint):
    '''Set of constraints to determine the value of DeviationVariables'''
    
    prefix = 'DVC_'
    
class HelperConstraint1(EnzymeConstraint):
    '''Set of constraints defined to change absolute value to linear terms in 
    objective function.'''
    
    prefix = 'HC1_'
    
class HelperConstraint2(EnzymeConstraint):
    '''Set of constraints defined to change absolute value to linear terms in 
    objective function.'''
    
    prefix = 'HC2_'


def mopa(perturb_model, solution, linear=False, keep_model=True, scaled=True,
         normalized=False):
    '''Compute a single solution based on (linear) MOPA (equivalent MOMA for ETFL).

    Compute a new flux distribution while trying to minimize the change of
    expression for the enzymes. Minimization of protein adjustment (MOPA)
    used to assess the impact of knock-outs. Thus the typical usage is to 
    provide a wildtype flux distribution as reference and a model in 
    knock-out state.

    Parameters
    ----------
    model : etfl.Model
        The model state to compute a MOMA-based solution for.
    solution : etfl.Solution, non-optional
        A (wildtype) reference solution.
    linear : bool, optional
        Whether to use the linear MOMA formulation or not (default True).
    keep_model_bkup : bool, optional
        Whether to keep the current model_bkup or change it for MOMA formulation
    scaled : bool, optional
        Whether to used real concentrations or the scaled ones (divided by molecular weight)
    normalized : bool, optional
        Whether to normalize deviations by reference values.

    Returns
    -------
    etfl.Solution
        A flux distribution that is at a minimal distance compared to the
        reference solution.
        
     Notes
    -----
    In the original MOMA [1]_ specification one looks for the flux distribution
    of the deletion (v^d) closest to the fluxes without the deletion (v). But,
    here we try to minimize the difference in enzyme concentration (E^d & E). 
    In math this means:

    minimize \sum_j (E^d_i - E_i)^2
    s.t. ETFL cons.


    In the case of linear MOPA, we instead minimize \sum_i abs(E^d_i - E_i). 
    The linear MOMA is typically significantly faster. Also quadratic MOMA tends
    to give flux distributions in which all fluxes deviate from the reference
    fluxes a little bit whereas linear MOMA tends to give flux distributions
    where the majority of fluxes are the same reference with few fluxes
    deviating a lot (typical effect of L2 norm vs L1 norm).

    '''
     
    if keep_model:
        model_bkup = deepcopy(perturb_model)
    else:
        model_bkup = perturb_model
        
    
    if not linear:
        model_bkup.solver = sutil.choose_solver(model_bkup, qp=True)
        
    # keeping the wild type growth rate
    setattr(model_bkup, 'original_growth', solution.objective_value)
    
    # resetting the objective function
    model_bkup.objective = model_bkup.problem.Objective(Zero, direction="min", sloppy=True)
    
    
    obj_vars = []
    str_modifier='EZ_{}'
    for enz in model_bkup.enzymes:
        try:
            this_conc = solution.raw.loc[str_modifier.format(enz.id)]
        except KeyError: # a protein migh be added to the model_perturbed which is not in the wild-type
            this_conc = 0
        
        coef = 1
        if not scaled:
            this_conc = this_conc / enz.scaling_factor
            coef = 1 / enz.scaling_factor
        this_var = model_bkup.add_variable(kind = DeviationVariable, 
                                      hook = enz, 
                                      lb = -1000, 
                                      ub = 1000)
        
        if normalized:
            if np.isclose(this_conc,0):
                sc_fact = 1/1e-3    # Basal flux
            else:
                sc_fact = 1/this_conc
            model_bkup.add_constraint(kind = DeviationConstraint, 
                             hook = enz, 
                             expr = this_var + \
                             model_bkup.variables[str_modifier.format(enz.id)]\
                             *coef*sc_fact ,
                             lb = 1, 
                             ub = 1)
        else:
            model_bkup.add_constraint(kind = DeviationConstraint, 
                             hook = enz, 
                             expr = this_var + \
                             model_bkup.variables[str_modifier.format(enz.id)]*coef ,
                             lb = this_conc, 
                             ub = this_conc)

        if linear:
            # In this case we need to add one set of variables and two sets of 
            # constraints.
            helper_var = model_bkup.add_variable(kind = HelperVariable, 
                                                 hook = enz, 
                                                 lb = 0, 
                                                 ub = 1000)
            model_bkup.add_constraint(kind = HelperConstraint1, 
                                      hook = enz, 
                                      expr = this_var - helper_var,
                                      lb = -1000, 
                                      ub = 0)
            model_bkup.add_constraint(kind = HelperConstraint2, 
                                      hook = enz, 
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
    mopa_sol = model_bkup.optimize()
    
        
    return mopa_sol