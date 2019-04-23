# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

ME-related Reaction subclasses and methods definition


"""

import pandas as pd
import numpy as np
from tqdm import tqdm

from ..optim.variables import EnzymeRef, mRNARef, EnzymeVariable, mRNAVariable, \
    LinearizationVariable
from ..optim.constraints import EnzymeDeltaPos, EnzymeDeltaNeg, \
    mRNADeltaPos, mRNADeltaNeg
from ..optim.utils import fix_growth, release_growth, \
                            get_active_growth_bounds, safe_optim

import operator

from pytfa.analysis.chebyshev import chebyshev_center

from collections import defaultdict


mrna_length_avg = 370
DEFAULT_DYNAMIC_CONS = {
    'mRNA_degradation':True,
    'mRNA_synthesis':True,
    'enzyme_synthesis':True,
    'enzyme_degradation':True,
}


class EnzymeDeltaRHS(EnzymeVariable):
    """
    E(t+dt) -  E(t) <= f(dt, E(t), mu)
    E(t+dt) <= E(t) + f(dt, E(t), mu)
    E(t+dt) <= ERHS(t, E(t), mu)
    """

    prefix = 'ERHS_'

class mRNADeltaRHS(mRNAVariable):
    """
    F(t+dt) -  F(t) <= f(dt, F(t), mu)
    F(t+dt) <= F(t) + f(dt, F(t), mu)
    F(t+dt) <= FRHS(t, F(t), mu)
    """

    prefix='FRHS_'

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

def add_enzyme_rhs_variable(dmodel):
    for enz in dmodel.enzymes:
        dmodel.add_variable(kind = EnzymeDeltaRHS,
                                hook = enz,
                                lb = 0,
                                ub = enz.variable.ub,
                                queue=True)

    dmodel.repair()

def add_mRNA_rhs_variable(dmodel):
    for mrna in dmodel.mrnas:
        dmodel.add_variable(kind = mRNADeltaRHS,
                                hook = mrna,
                                lb = 0,
                                ub = mrna.variable.ub,
                                queue=True)

    dmodel.repair()
    
def get_mu_times_var(dmodel, macromolecule):
    resolution = dmodel.mu_bins[-1][-1][-1] / len(dmodel.mu_bins)

    growth_discretization_vars = dmodel.get_ordered_ga_vars()

    # Build z =   ga_0*2^0*mu_max/N * [E]
    #           + ga_1*2^1*mu_max/N * [E]
    #           + ...
    #           + ga_n*2^n*mu_max/N * [E]
    out_expr = 0
    for i, ga_i in enumerate(growth_discretization_vars):
        # Linearization step for ga_i * [E]
        z_name = LinearizationVariable.prefix + \
                 '__MUL__'.join([ga_i.name,macromolecule.variable.name])
        model_z_i = dmodel.variables.get(z_name)
        out_expr += (2 ** i) * resolution * model_z_i

    return out_expr

def add_enzyme_delta_constraint(dmodel, timestep, degradation, synthesis):
    """
    Adds the constraint

    E-Eref <= Δt*v_assembly_max
    E-Eref-Δt*v_assembly_max <= 0


    :param dmodel:
    :param timestep:
    :return:
    """

    enzyme_ref_variables = dmodel.get_variables_of_type(EnzymeRef)
    # enzyme_rhs_variables = dmodel.get_variables_of_type(EnzymeDeltaRHS)

    for enz in dmodel.enzymes:

        if 'dummy' in enz.id:
            continue

        # E = enz.concentration
        # E_ref = enzyme_ref_variables.get_by_id(enz.id) * enz.scaling_factor
        # E_rhs = enzyme_rhs_variables.get_by_id(enz.id)
        E_ref = enzyme_ref_variables.get_by_id(enz.id)

        v_asm = (enz.complexation.forward_variable \
                    - enz.complexation.reverse_variable) \
                * enz.complexation.scaling_factor \
                / enz.scaling_factor

        # Backwards Euler method
        # E(t+dt) = (mu(t+dt)+kdeg)*E(t+dt)*dt + E(t)
        E = enz.variable
        # muE = get_mu_times_var(dmodel, enz)
        # conc_change = (enz.kdeg * E + muE) * timestep
        conc_change = v_asm * timestep

        if synthesis:
            # expr_pos = (E - E_ref) - timestep*v_asm
            # expr_pos = E - E_rhs
            expr_pos = E - E_ref - conc_change

            dmodel.add_constraint(kind = EnzymeDeltaPos,
                                      hook = enz,
                                      expr = expr_pos,
                                      ub = 0,
                                      queue = True
                                      )
        if degradation:
            # expr_neg = (E_ref - E) - timestep*v_asm
            # expr_neg = E_rhs - E
            expr_neg = E_ref - E - conc_change

            dmodel.add_constraint(kind = EnzymeDeltaNeg,
                                      hook = enz,
                                      expr = expr_neg,
                                      ub = 0,
                                      queue = True
                                      )
    dmodel.repair()


def add_mRNA_delta_constraint(dmodel, timestep, degradation, synthesis):
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
    # mrna_rhs_variables = dmodel.get_variables_of_type(mRNADeltaRHS)


    for mrna in dmodel.mrnas:

        if 'dummy' in mrna.id:
            continue

        # F = mrna.concentration
        F_ref = mrna_ref_variables.get_by_id(mrna.id)
        # F_rhs = mrna_rhs_variables.get_by_id(mrna.id)

        trans_id = dmodel._get_transcription_name(mrna.id)

        trans = dmodel.transcription_reactions.get_by_id(trans_id)
        v_syn = (trans.forward_variable - trans.reverse_variable) \
                 * trans.scaling_factor / mrna.scaling_factor

        # Backwards Euler method
        # F(t+dt) = (mu(t+dt)+kdeg)*F(t+dt)*dt + F(t)
        F = mrna.variable
        # muF = get_mu_times_var(dmodel, mrna)
        # conc_change = (mrna.kdeg * F + muF) * timestep
        conc_change = v_syn * timestep


        if synthesis:
            # expr_pos = (F - F_ref)- timestep*v_syn
            # expr_pos = F - F_rhs
            expr_pos = F - F_ref - conc_change

            dmodel.add_constraint(kind = mRNADeltaPos,
                                      hook = mrna,
                                      expr = expr_pos,
                                      ub = 0,
                                      queue = True
                                      )

        if degradation:
            # expr_neg = (F_ref - F) - timestep*v_syn
            # expr_neg = F_rhs - F
            expr_neg = F_ref - F - conc_change

            dmodel.add_constraint(kind = mRNADeltaNeg,
                                      hook = mrna,
                                      expr = expr_neg,
                                      ub = 0,
                                      queue = True
                                      )
    dmodel.repair()



def add_dynamic_variables_constraints(dmodel, timestep, dynamic_constraints):
    if dynamic_constraints['mRNA_degradation'] or dynamic_constraints['mRNA_synthesis']:
        add_mRNA_ref_variable(dmodel)
        # add_mRNA_rhs_variable(dmodel)
        add_mRNA_delta_constraint(dmodel, timestep,
                                  degradation=dynamic_constraints['mRNA_degradation'],
                                  synthesis  =dynamic_constraints['mRNA_synthesis'],
                                  )
    if dynamic_constraints['enzyme_degradation'] or dynamic_constraints['enzyme_synthesis']:
        add_enzyme_ref_variable(dmodel)
        # add_enzyme_rhs_variable(dmodel)
        add_enzyme_delta_constraint(dmodel, timestep,
                                  degradation=dynamic_constraints['enzyme_degradation'],
                                  synthesis  =dynamic_constraints['enzyme_synthesis'],
                                  )


def apply_ref_state(dmodel, solution, timestep, has_mrna, has_enzymes):

    enz_ref  = dmodel.get_variables_of_type(EnzymeRef)
    mrna_ref = dmodel.get_variables_of_type(mRNARef)
    # enz_rhs  = dmodel.get_variables_of_type(EnzymeDeltaRHS)
    # mrna_rhs = dmodel.get_variables_of_type(mRNADeltaRHS)

    epsilon = dmodel.solver.configuration.tolerances.feasibility

    if has_enzymes:
        for enz in dmodel.enzymes:
            # the_rhs = enz_rhs.get_by_id(enz.id)
            the_ref = enz_ref.get_by_id(enz.id)
            E0 = solution.loc[enz.variable.name]
            the_value = E0
            # Taylor expansion element
            # e = (enz.kdeg + solution.loc[dmodel.growth_reaction.id]) * timestep
            # the_value = E0 * (1 + e/1 + (e**2)/2 + (e**3)/6)

            try:
                # the_rhs.variable.ub = max(0,the_value + epsilon)
                # the_rhs.variable.lb = max(0,the_value - epsilon)
                the_ref.variable.ub = max(0,the_value + epsilon)
                the_ref.variable.lb = max(0,the_value - epsilon)
            except ValueError:
                # the_rhs.variable.lb = max(0,the_value - epsilon)
                # the_rhs.variable.ub = max(0,the_value + epsilon)
                the_ref.variable.lb = max(0,the_value - epsilon)
                the_ref.variable.ub = max(0,the_value + epsilon)
    if has_mrna:
        for mrna in dmodel.mrnas:
            # the_rhs = mrna_rhs.get_by_id(mrna.id)
            the_ref = mrna_ref.get_by_id(mrna.id)
            F0 = solution.loc[mrna.variable.name]
            # Taylor expansion element
            # e = (mrna.kdeg + solution.loc[dmodel.growth_reaction.id]) * timestep
            # the_value = F0 * (1 + e/1 + (e**2)/2 + (e**3)/6)
            the_value = F0

            try:
                # the_rhs.variable.ub = max(0,the_value + epsilon)
                # the_rhs.variable.lb = max(0,the_value - epsilon)
                the_ref.variable.ub = max(0,the_value + epsilon)
                the_ref.variable.lb = max(0,the_value - epsilon)
            except ValueError:
                # the_rhs.variable.lb = max(0,the_value - epsilon)
                # the_rhs.variable.ub = max(0,the_value + epsilon)
                the_ref.variable.lb = max(0,the_value - epsilon)
                the_ref.variable.ub = max(0,the_value + epsilon)

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
    # dmodel.growth_reaction.lower_bound = mu_lb - 0.0001

    fix_growth(dmodel, dmodel.solution)

    try:
        dmodel.objective = dmodel.chebyshev_radius.radius.variable
        dmodel.optimize()
        chebyshev_sol = dmodel.solution

        release_growth(dmodel)
    except AttributeError: #does not solve
        dmodel.logger.warning('Chebyshev solution at fixed growth solution '
                              'not found - relaxing integers and solving '
                              'with growth lb')
        # We relax the integer bounds and try solving anyway with just the lb
        epsilon = dmodel.solver.configuration.tolerances.feasibility
        release_growth(dmodel)
        dmodel.growth_reaction.lower_bound = mu_lb - epsilon
        dmodel.optimize()

    dmodel.growth_reaction.lower_bound = prev_lb
    dmodel.objective = prev_obj
    return chebyshev_sol


def show_initial_solution(model, solution):
    print('Objective            : {}'.format(solution.objective_value))
    print(' - Growth            : {}'.format(solution.raw.loc[model.growth_reaction.id]))
    print(' - Ribosomes produced: {}'.format(solution.raw.loc[model.ribosome.variable.name]))
    print(' - RNAP produced: {}'.format(solution.raw.loc[model.rnap.variable.name]))


BIGM=1000

def run_dynamic_etfl(model, timestep, tfinal, uptake_fun, medium_fun,
                     uptake_enz,
                     S0, X0, inplace=False, initial_solution = None,
                     chebyshev_bigm=BIGM, chebyshev_variables=None,
                     chebyshev_exclude=None, chebyshev_include=None,
                     dynamic_constraints=DEFAULT_DYNAMIC_CONS):
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
        chebyshev_exclude = []#\
        # [LinearizationVariable,
        #                  EnzymeDeltaNeg,
        #                  EnzymeDeltaPos,
        #                  mRNADeltaNeg,
        #                  mRNADeltaPos]

    if chebyshev_include is None:
        chebyshev_include = list()

    if chebyshev_variables is None:
        from ..optim.variables import mRNAVariable, EnzymeVariable
        chebyshev_variables =  dmodel.get_variables_of_type(mRNAVariable)
        chebyshev_variables += dmodel.get_variables_of_type(EnzymeVariable)


    add_dynamic_variables_constraints(dmodel, timestep, dynamic_constraints)
    chebyshev_center(dmodel, chebyshev_variables,
                     inplace=True,
                     big_m=chebyshev_bigm,
                     include=chebyshev_include,
                     exclude=chebyshev_exclude)
    dmodel.objective = the_obj

    chebyshev_sol = compute_center(dmodel)


    has_mrna = dynamic_constraints['mRNA_degradation'] + dynamic_constraints['mRNA_synthesis']
    has_enzymes = dynamic_constraints['enzyme_degradation'] + dynamic_constraints['enzyme_synthesis']

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
    dmodel.optimize()

    show_initial_solution(dmodel, current_solution)

    for  k, t in tqdm(enumerate(times)):

        apply_ref_state(dmodel, current_solution.raw, timestep, has_mrna, has_enzymes)

        for uptake_flux, kinfun in uptake_fun.items():
            # Concentration C, cell density X
            # C(t+dt) = C(t) - dt*v*X
            # C(t+dt) >= 0 => v <= C(t)/(dt*X)
            max_available_substrate = S[uptake_flux] / (timestep * X)

            if uptake_flux in uptake_enz:
                these_uptake_enz = [enz for x in uptake_enz[uptake_flux] for enz in dmodel.reactions.get_by_id(x).enzymes]
                vmax = sum([x.kcat_fwd * x.scaling_factor * dmodel.solution.raw[x.variable.name] for x in these_uptake_enz])

                lb = vmax * kinfun(S[uptake_flux])
                dmodel.reactions.get_by_id(uptake_flux).lower_bound = -1* min(lb, max_available_substrate)

            else:
                lb = kinfun(S[uptake_flux])
                dmodel.reactions.get_by_id(uptake_flux).lower_bound = -1* min(lb, max_available_substrate)



        try:
            the_solution = compute_center(dmodel)
        except AttributeError:
            print('############################')
            print('### Crashed at t={},k={}'.format(t,k))
            print('############################')
            return wrap_time_sol(var_solutions, obs_values)

        colname = 't_{}'.format(k)

        var_solutions[colname] = the_solution.raw.copy()
        X,S= update_medium(t,X,S,model,medium_fun,timestep)
        obs_values = update_sol(t,X,S,dmodel,obs_values, colname)

        current_solution = the_solution
        wrap_time_sol(var_solutions, obs_values).to_csv('tmp_detfl.csv')

    return wrap_time_sol(var_solutions, obs_values)

def wrap_time_sol(var_solutions, obs_values):
    return pd.concat([var_solutions,obs_values], axis=0)


