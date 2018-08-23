# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: ETFL team

ME-related Reaction subclasses and methods definition


"""

import pandas as pd
from ..optim.variables import EnzymeRef, mRNARef
from ..optim.constraints import EnzymeDeltaPos, EnzymeDeltaNeg, \
    mRNADeltaPos, mRNADeltaNeg
import operator

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
        E = enz.variable
        E_ref = enzyme_ref_variables.get_by_id(enz.id)

        v_loss = (enz.complexation.forward_variable \
                    - enz.complexation.reverse_variable)

        # Find length(pep)*stoichiometry for all peptides

        peptides = enz.complexation.metabolites
        pep_weighted_len =  {
            pep:len(dmodel.peptides.get_by_id(pep.id).peptide)*(-1)*stoich
            for pep,stoich in peptides.items()
        }

        maxpep, coeff = max(pep_weighted_len.items(), key=operator.itemgetter(1))
        rib_usage = dmodel.ribosome_usage.get_by_id(maxpep.id)

        # v_asm = dmodel.ribosome.kribo \
        #         * 1/coeff \
        #         * rib_usage

        v_asm = dmodel.ribosome.kribo \
                * 1/sum(pep_weighted_len.values()) \
                * (dmodel.Rt - dmodel.Rf)

        # expr_pos = (E - E_ref) - timestep*v_asm
        expr_pos = (E - E_ref) - timestep*v_loss

        dmodel.add_constraint(kind = EnzymeDeltaPos,
                                  hook = enz,
                                  expr = expr_pos,
                                  ub = 0,
                                  queue = True
                                  )

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
        F = mrna.variable
        F_ref = mrna_ref_variables.get_by_id(mrna.id)
        trans_id = dmodel._get_transcription_name(mrna.id)
        if trans_id == 'dummy_gene_transcription':
            trans_id = 'dummy_transcription'
        trans = dmodel.transcription_reactions.get_by_id(trans_id)
        v_loss = (trans.forward_variable - trans.reverse_variable)

        mrna_length = len(mrna.rna)
        v_syn = dmodel.rnap.ktrans*dmodel.rnap.variable/mrna_length \
                                      *dmodel._mrna_scaling/dmodel._prot_scaling

        # expr_pos = (F - F_ref)- timestep*v_syn
        expr_pos = (F - F_ref)- timestep*v_loss

        dmodel.add_constraint(kind = mRNADeltaPos,
                                  hook = mrna,
                                  expr = expr_pos,
                                  ub = 0,
                                  queue = True
                                  )

        expr_neg = (F_ref - F) - timestep*v_loss
        # expr_neg = (F_ref - F) - timestep*v_syn

        dmodel.add_constraint(kind = mRNADeltaNeg,
                                  hook = mrna,
                                  expr = expr_neg,
                                  ub = 0,
                                  queue = True
                                  )
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


def run_dynamic_etfl(model, timestep, tfinal, uptake_fun, medium_fun,
                     S0, X0, inplace=False):
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

    add_dynamic_variables_constraints(dmodel, timestep)

    # add_dynamic_boundaries(dmodel)

    t = 0.0

    S = S0
    X = X0

    St = defaultdict(list)
    upt_lim = defaultdict(list)
    upt = defaultdict(list)
    Xt = []
    tt = []
    mut = []


    for uptake_flux, kinfun in uptake_fun.items():
        the_rxn = dmodel.reactions.get_by_id(uptake_flux)
        the_rxn.upper_bound = 0
        the_rxn.lower_bound = -1 * kinfun(S0[uptake_flux])


    # dmodel.reactions.EX_o2_e.lower_bound = -15 # Mahadevan et al. 2002
    dmodel.optimize()

    time_solution = pd.DataFrame()

    k = 0
    while t<tfinal:
        apply_ref_state(dmodel, dmodel.solution.x_dict)

        for uptake_flux, kinfun in uptake_fun.items():
            dmodel.reactions.get_by_id(uptake_flux).lower_bound = \
                                        -1 * kinfun(S[uptake_flux])

        try:
            dmodel.optimize()
        except AttributeError:
            # raise(AttributeError)
            time_solution = wrap_time_sol(time_solution, tt, Xt, St, mut, upt, upt_lim)
            return time_solution

        mu = dmodel.growth_reaction.flux

        time_solution['t_{}'.format(k)] = dmodel.solution.x_dict

        t += timestep
        Xt.append(X)
        tt.append(t)
        mut.append(mu)
        X += mu * X * timestep

        for uptake_flux, medium_change in medium_fun.items():
            sol_flux = dmodel.reactions.get_by_id(uptake_flux).flux
            St[uptake_flux].append(S[uptake_flux])
            upt_lim[uptake_flux].append(uptake_fun[uptake_flux](S[uptake_flux]))
            upt[uptake_flux].append(-1*sol_flux)
            S[uptake_flux] += sol_flux * X * timestep
            S[uptake_flux] = medium_change(t,S[uptake_flux])

        k+=1

    time_solution = wrap_time_sol(time_solution, tt, Xt, St, mut, upt, upt_lim)
    return time_solution


def wrap_time_sol(time_solution, tt, Xt, St, mut, upt, upt_lim):

    for uptake_flux, values in St.items():
        time_solution.loc['S_'+uptake_flux] = values
        time_solution.loc['ub_'+uptake_flux] = upt[uptake_flux]
        time_solution.loc['act_'+uptake_flux] = upt_lim[uptake_flux]
    time_solution.loc['X'] = Xt
    time_solution.loc['t'] = tt
    time_solution.loc['mu'] = mut

    return time_solution


