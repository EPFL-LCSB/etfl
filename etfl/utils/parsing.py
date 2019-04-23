"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

Parsing utilities

"""

import re
import sympy
from sympy.parsing.sympy_parser import parse_expr
import ast

ESCAPE_CHARS = (
    r'\n',
    r'\t',
    r'\s',
    r'"',
    r"'",
)

def isevaluable(s):
    """
    Test evaluability of a string for eval with sympy

    :param s:
    :return:
    """
    # try:
    #     ast.literal_eval(s)
    #     return True
    # except ValueError:
    #     return False

    return not any(x in s for x in ESCAPE_CHARS)

def parse_gpr(gpr):
    """
    Parses a string gpr into a sympy expression

    :param gpr:
    :return:
    """

    #assert(isinstance(gpr,str))
    assert(isevaluable(gpr))

    #TODO: Faster implementation with single pass regex ?
    bool_gpr = gpr.lower().replace('or','|').replace('and','&')

    sym_gpr = parse_expr(bool_gpr)

    return sym_gpr