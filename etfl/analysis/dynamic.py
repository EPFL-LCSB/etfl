# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: ETFL team

ME-related Reaction subclasses and methods definition


"""

import pandas as pd
import numpy as np
from tqdm import tqdm

from ..optim.variables import EnzymeRef, mRNARef
from ..optim.constraints import EnzymeDeltaPos, EnzymeDeltaNeg, \
    mRNADeltaPos, mRNADeltaNeg
import operator

from pytfa.analysis.chebyshev import chebyshev_center

from collections import defaultdict

mrna_length_avg = 370

def add_enzyme_ref_variable(dmodel):
    for enz in dmodel.enzymes:
        dmodel.add_variable(kind = EnzymeRef,
                                hook = enz,
                                lb = 0,
                                ub = enz.variable.ub,
                                queue=True)

    dmodel.repair()

def add_mRNA_ref_variable(dmodel):
    for mrna in dmodel.mrnas:
        dmodel.add_variable(kind = mRNARef,
                                hook = mrna,
                                lb = 0,
                                ub = mrna.variable.ub,
                                queue=True)

    dmodel.repair()

def add_enzyme_delta_constraint(dmodel, timestep):
    """
    Adds the constraint

    E-Eref <= Δt*v_assembly_max
    E-Eref-Δt*v_assembly_max <= 0


    :param dmodel:
    :param timestep:
    :return:
    """

    enzyme_ref_variables = dmodel.get_variables_of_type(EnzymeRef)

    for enz in dmodel.enzymes:

        if 'dummy' in enz.id:
            continue

        E = enz.concentration
        E_ref = enzyme_ref_variables.get_by_id(enz.id) * enz.scaling_factor

        v_loss = (enz.complexation.forward_variable \
                    - enz.complexation.reverse_variable) \
                 * enz.complexation.scaling_factor

        # expr_pos = (E - E_ref) - timestep*v_asm
        expr_pos = (E - E_ref) - timestep*v_loss

        # dmodel.add_constraint(kind = EnzymeDeltaPos,
        #                           hook = enz,
        #                           expr = expr_pos,
        #                           ub = 0,
        #                           queue = True
        #                           )

        expr_neg = (E_ref - E) - timestep*v_loss
        # expr_neg = (E_ref - E) - timestep*v_asm

        dmodel.add_constraint(kind = EnzymeDeltaNeg,
                                  hook = enz,
                                  expr = expr_neg,
                                  ub = 0,
                                  queue = True
                                  )
    dmodel.repair()


def add_mRNA_delta_constraint(dmodel, timestep):
    """
    Adds the constraint

        |F-Fref| <= Δt*v_transcription_max
    <=>
        F-Fref -Δt*v_transcription_max <= 0
        &
        Fref-F -Δt*v_transcription_max <= 0



    :param dmodel:
    :param timestep:
    :return:
    """

    mrna_ref_variables = dmodel.get_variables_of_type(mRNARef)

    for mrna in dmodel.mrnas:

        if 'dummy' in mrna.id:
            continue

        F = mrna.concentration
        F_ref = mrna_ref_variables.get_by_id(mrna.id) * mrna.scaling_factor
        trans_id = dmodel._get_transcription_name(mrna.id)

        trans = dmodel.transcription_reactions.get_by_id(trans_id)
        v_loss = (trans.forward_variable - trans.reverse_variable) \
                 * trans.scaling_factor


        expr_pos = (F - F_ref)- timestep*v_loss

        # dmodel.add_constraint(kind = mRNADeltaPos,
        #                           hook = mrna,
        #                           expr = expr_pos,
        #                           ub = 0,
        #                           queue = True
        #                           )
        #
        # expr_neg = (F_ref - F) - timestep*v_loss
        # # expr_neg = (F_ref - F) - timestep*v_syn
        #
        # dmodel.add_constraint(kind = mRNADeltaNeg,
        #                           hook = mrna,
        #                           expr = expr_neg,
        #                           ub = 0,
        #                           queue = True
        #                           )
    dmodel.repair()



def add_dynamic_variables_constraints(dmodel, timestep):
    add_enzyme_ref_variable(dmodel)
    add_enzyme_delta_constraint(dmodel, timestep)
    add_mRNA_ref_variable(dmodel)
    add_mRNA_delta_constraint(dmodel, timestep)

def add_dynamic_boundaries(dmodel):
    dmodel.reactions.EX_glc__D_e.lower_bound = -10

def apply_ref_state(dmodel, solution):

    enz_ref  = dmodel.get_variables_of_type(EnzymeRef)
    mrna_ref = dmodel.get_variables_of_type(mRNARef)

    epsilon = dmodel.solver.configuration.tolerances.feasibility

    for enz in dmodel.enzymes:
        the_ref = enz_ref.get_by_id(enz.id)

        try:
            the_ref.variable.ub = max(0,solution[enz.variable.name] + epsilon)
            the_ref.variable.lb = max(0,solution[enz.variable.name] - epsilon)
        except ValueError:
            the_ref.variable.lb = max(0,solution[enz.variable.name] - epsilon)
            the_ref.variable.ub = max(0,solution[enz.variable.name] + epsilon)

    for mrna in dmodel.mrnas:
        the_ref = mrna_ref.get_by_id(mrna.id)

        try:
            the_ref.variable.ub = max(0,solution[mrna.variable.name] + epsilon)
            the_ref.variable.lb = max(0,solution[mrna.variable.name] - epsilon)
        except ValueError:
            the_ref.variable.lb = max(0,solution[mrna.variable.name] - epsilon)
            the_ref.variable.ub = max(0,solution[mrna.variable.name] + epsilon)

def update_sol(t, X, S, dmodel, obs_values, colname):
    obs_values.loc['t',colname] = t
    obs_values.loc['X',colname] = X
    obs_values.loc['mu',colname] = dmodel.solution.fluxes[dmodel.growth_reaction.id]

    obs_values.loc['X', colname] = X
    for uptake_flux, value in S.items():
        the_rxn = dmodel.reactions.get_by_id(uptake_flux)
        obs_values.loc['S_'+uptake_flux, colname] = value
        obs_values.loc['ub_'+uptake_flux, colname] = -1*the_rxn.lower_bound
        obs_values.loc['act_'+uptake_flux, colname] = -1*dmodel.solution.fluxes[uptake_flux]
    return obs_values

def update_medium(t, Xi, Si, dmodel, medium_fun, timestep):
    mu = dmodel.solution.fluxes[dmodel.growth_reaction.id]
    X = Xi + mu * Xi * timestep
    S = dict()
    for uptake_flux, medium_change in medium_fun.items():
        sol_flux = dmodel.solution.fluxes[uptake_flux]
        S[uptake_flux] = Si[uptake_flux] + sol_flux * Xi * timestep
        S[uptake_flux] = medium_change(t, S[uptake_flux])

    return X,S

def get_active_growth_bounds(model):
    mu = model.growth_reaction.flux
    difflist = [abs(mu - x[0]) for x in model.mu_bins]
    min_diff = min(difflist)
    min_ix = difflist.index(min_diff)

    mu_i, (mu_lb, mu_ub) = model.mu_bins[min_ix]

    return mu_i, mu_lb, mu_ub

BIGM=1000


def compute_center(dmodel):
    """
    Fixes growth to be above computed lower bound, finds chebyshev center,
    resets the model, returns solution data

    :param dmodel:
    :return:
    """
    dmodel.optimize() # Almost free if the model has been optimized before
    prev_lb = dmodel.growth_reaction.lower_bound
    prev_obj = dmodel.objective.expression
    _,mu_lb,_ = get_active_growth_bounds(dmodel)
    dmodel.growth_reaction.lower_bound = mu_lb - 0.0001
    dmodel.objective = dmodel.chebyshev_radius.radius.variable
    dmodel.optimize()
    chebyshev_sol = dmodel.solution
    dmodel.growth_reaction.lower_bound = prev_lb
    dmodel.objective = prev_obj
    return chebyshev_sol



def run_dynamic_etfl(model, timestep, tfinal, uptake_fun, medium_fun,
                     S0, X0, inplace=False, initial_solution = None,
                     chebyshev_bigm=BIGM, chebyshev_variables=None,
                     chebyshev_exclude=None, chebyshev_include=None):
    """

    :param model:
    :param initial_solution:
    :param timestep:
    :param tfinal:
    :return:
    """

    if not inplace:
        dmodel = model.copy()
    else:
        dmodel = model

    dmodel.optimize()

    the_obj = dmodel.objective.expression

    # Add chebyshev center transformation to stay in an approximate center of
    # the feasible space at all steps

    if chebyshev_exclude is None:
        from ..optim.variables import LinearizationVariable
        chebyshev_exclude = [LinearizationVariable]

    if chebyshev_include is None:
        chebyshev_include = list()

    if chebyshev_variables is None:
        from ..optim.variables import mRNAVariable, EnzymeVariable
        chebyshev_variables =  dmodel.get_variables_of_type(mRNAVariable)
        chebyshev_variables += dmodel.get_variables_of_type(EnzymeVariable)

    chebyshev_center(dmodel, chebyshev_variables,
                     inplace=True,
                     big_m=chebyshev_bigm,
                     include=chebyshev_include,
                     exclude=chebyshev_exclude)
    dmodel.objective = the_obj

    chebyshev_sol = compute_center(dmodel)

    add_dynamic_variables_constraints(dmodel, timestep)

    for uptake_flux, kinfun in uptake_fun.items():
        the_rxn = dmodel.reactions.get_by_id(uptake_flux)
        # the_rxn.upper_bound = 0
        the_rxn.lower_bound = -1 * kinfun(S0[uptake_flux])

    if initial_solution is None:
        current_solution = chebyshev_sol
    else:
        current_solution = initial_solution

    var_solutions = pd.DataFrame()
    obs_values = pd.DataFrame()

    times = np.arange(0,tfinal,timestep)
    S = S0
    X = X0

    for  k, t in tqdm(enumerate(times)):

        apply_ref_state(dmodel, current_solution.raw)

        for uptake_flux, kinfun in uptake_fun.items():
            dmodel.reactions.get_by_id(uptake_flux).lower_bound = \
                -1 * kinfun(S[uptake_flux])

        try:
            the_solution = compute_center(dmodel)
        except AttributeError:
            print('############################')
            print('### Crashed at t={},k={}'.format(t,k))
            print('############################')
            return wrap_time_sol(var_solutions, obs_values)

        colname = 't_{}'.format(k)

        var_solutions[colname] = the_solution.raw.copy()
        obs_values = update_sol(t,X,S,dmodel,obs_values, colname)
        X,S= update_medium(t,X,S,model,medium_fun,timestep)

        current_solution = the_solution

    return wrap_time_sol(var_solutions, obs_values)

def wrap_time_sol(var_solutions, obs_values):
    return pd.concat([var_solutions,obs_values], axis=0)


