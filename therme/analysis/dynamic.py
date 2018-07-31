# -*- coding: utf-8 -*-
"""
.. module:: thermome
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: ETFL team

ME-related Reaction subclasses and methods definition


"""

import pandas as pd
from ..optim.variables import EnzymeRef, mRNARef
from ..optim.constraints import EnzymeDeltaPos, EnzymeDeltaNeg, \
    mRNADeltaPos, mRNADeltaNeg

def add_enzyme_ref_variable(dmodel):
    for enz in dmodel.enzymes:
        dmodel.add_variable(kind = EnzymeRef,
                                hook = enz,
                                lb = enz.variable.lb,
                                ub = enz.variable.ub,
                                queue=False)

    dmodel._update()
    dmodel.repair()

def add_mRNA_ref_variable(dmodel):
    for mrna in dmodel.mrnas:
        dmodel.add_variable(kind = mRNARef,
                                hook = mrna,
                                lb = mrna.variable.lb,
                                ub = mrna.variable.ub,
                                queue=False)

    dmodel._update()
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
        v_asm = enz.complexation.forward_variable \
                - enz.complexation.reverse_variable

        expr_pos = (E - E_ref) - timestep*v_asm

        dmodel.add_constraint(kind = EnzymeDeltaPos,
                                  hook = enz,
                                  expr = expr_pos,
                                  ub = 0,
                                  queue = False
                                  )

        expr_neg = (E_ref - E) - timestep*v_asm

        dmodel.add_constraint(kind = EnzymeDeltaNeg,
                                  hook = enz,
                                  expr = expr_neg,
                                  ub = 0,
                                  queue = False
                                  )
    dmodel._update()
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
        trans = dmodel.transcription_reactions.get_by_id(trans_id)
        v_syn = trans.forward_variable - trans.reverse_variable

        expr_pos = (F - F_ref)- timestep*v_syn

        dmodel.add_constraint(kind = mRNADeltaPos,
                                  hook = mrna,
                                  expr = expr_pos,
                                  ub = 0,
                                  queue = True
                                  )

        expr_neg = (F_ref - F) - timestep*v_syn

        dmodel.add_constraint(kind = mRNADeltaNeg,
                                  hook = mrna,
                                  expr = expr_neg,
                                  ub = 0,
                                  queue = True
                                  )
    dmodel._update()
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
            the_ref.variable.ub = solution[enz.variable.name] + epsilon
            the_ref.variable.lb = solution[enz.variable.name] - epsilon
        except ValueError:
            the_ref.variable.lb = solution[enz.variable.name] - epsilon
            the_ref.variable.ub = solution[enz.variable.name] + epsilon

    for mrna in dmodel.mrnas:
        the_ref = mrna_ref.get_by_id(mrna.id)

        try:
            the_ref.variable.ub = solution[mrna.variable.name] + epsilon
            the_ref.variable.lb = solution[mrna.variable.name] - epsilon
        except ValueError:
            the_ref.variable.lb = solution[mrna.variable.name] - epsilon
            the_ref.variable.ub = solution[mrna.variable.name] + epsilon


def run_dynamic_etfl(model, initial_solution, timestep, tfinal):
    """

    :param model:
    :param initial_solution:
    :param timestep:
    :param tfinal:
    :return:
    """

    dmodel = model.copy()
    add_dynamic_variables_constraints(dmodel, timestep)

    add_dynamic_boundaries(dmodel)

    time_solution = pd.DataFrame(initial_solution)
    time_solution.columns = ['0']

    t = 0

    while t<tfinal:
        if t == 0:
            apply_ref_state(dmodel, initial_solution)
        else:
            apply_ref_state(dmodel, dmodel.solution.x_dict)

        dmodel.optimize()

        time_solution['{:.2}'.format(t+timestep)] = dmodel.solution.x_dict
        t+= timestep

    return time_solution


