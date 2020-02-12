# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

Solver configuration helpers

"""

def standard_solver_config(model, verbose = True):
    """
    Basic solver settings for ETFL
    :param model:
    :param verbose:
    :return:
    """
    # Solver settings
    model.solver.configuration.verbosity = int(verbose)
    model.solver.configuration.tolerances.feasibility = 1e-9
    model.solver.configuration.timeout = 7200
    # model.solver.configuration.presolve = True
    # model.solver.configuration.lp_method = 'dual'
    # model.solver.configuration.lp_method = 'barrier'
    if model.solver.interface.__name__ == 'optlang.gurobi_interface':
        # model.solver.problem.Params.BarHomogeneous = 1
        model.solver.problem.Params.NumericFocus = 3
        # model.solver.problem.Params.ScaleFlag = 3  # Aggressive
        model.solver.problem.Params.Method = -1  # Auto
        # model.solver.problem.Params.Method = 5  # deterministic concurrent simplex
        # model.solver.problem.Params.Aggregate = 0  # Turn off aggregation


def gene_ko_config(model):
    """
    Solver settings for performing gene KO. Tuned using the grbtune tool on the
    vETFL model iJO1366. The gene KO analysis is turned into a feasibility
    problem by putting a lower bound on growth.

    :param model:
    :return:
    """
    standard_solver_config(model)
    if model.solver.interface.__name__ == 'optlang.gurobi_interface':
        # model.solver.problem.Params.MIPFocus = 1# Find a feasible solution
        model.solver.problem.Params.NormAdjust = 2#
        model.solver.problem.Params.AggFill = 10#
        model.solver.problem.Params.BranchDir = 1#
        model.solver.problem.Params.PrePasses = 5#


def growth_uptake_config(model):
    """
    Solver settings for performing growth vs uptake studies. Tuned using the
    grbtune tool on the vETFL model iJO1366.

    :param model:
    :return:
    """
    standard_solver_config(model)
    if model.solver.interface.__name__ == 'optlang.gurobi_interface':
        # model.solver.problem.Params.Cuts = 3 # Focus on finding feasible solutions
        # model.solver.problem.Params.MIPFocus = 1 # Focus on finding feasible solutions
        # model.solver.problem.Params.ConcurrentMIP = 4 # 4 independant solvers
        model.solver.problem.Params.SimplexPricing = 2
        model.solver.problem.Params.BranchDir = 1
        model.solver.problem.Params.Heuristics = 0
        model.solver.problem.Params.AggFill = 5



def redhuman_config(model):
    """
    Solver settings for optimizing growth on human cancer models. Tuned using the
    grbtune tool on the vETFL model from reduced RECON3.

    :param model:
    :return:
    """
    standard_solver_config(model)
    if model.solver.interface.__name__ == 'optlang.gurobi_interface':
        model.solver.problem.Params.NormAdjust = 0
        model.solver.problem.Params.RINS = 100
        model.solver.problem.Params.ZeroObjNodes = 2500
        model.solver.problem.Params.NumericFocus = 3