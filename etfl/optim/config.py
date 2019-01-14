# -*- coding: utf-8 -*-
"""
.. module:: thermome
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: ETFL team

Solver configuration helpers

"""

def standard_solver_config(model, verbose = True):
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