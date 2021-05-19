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

GPR2EXPR_SUBS_DICT = {' or ': ' | ',
                      ' and ': ' & '}
EXPR2GPR_SUBS_DICT = {'|':' or ',
                      '&':' and '}

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

    bool_gpr = multiple_replace(gpr, GPR2EXPR_SUBS_DICT, ignore_case=True)
    sym_gpr = sympy.to_dnf(bool_gpr)

    return sym_gpr

def multiple_replace(text, adict, ignore_case = False):
    """ From https://www.oreilly.com/library/view/python-cookbook/0596001673/ch03s15.html
    """
    # the simplest, lambda-based implementation
    # Create a regular expression from all of the dictionary keys
    if ignore_case:
        regex = re.compile(r'\b' + r'\b|\b'.join(map(re.escape, adict.keys(  ))) + r'\b', re.IGNORECASE)
    else:
        regex = re.compile(r'\b' + r'\b|\b'.join(map(re.escape, adict.keys(  ))) + r'\b')

# For each match, look up the corresponding value in the dictionary
    return regex.sub(lambda match: adict[match.group(0)], text)


def simplify_gpr(gpr):
    from sympy import simplify
    formatted_gpr = gpr2expr(gpr)
    simplified_formatted_gpr = str(simplify(formatted_gpr))
    new_formatted_gpr = expr2gpr(simplified_formatted_gpr)

    print('Reduced expression from {} char. to {}: {}'.format(len(gpr),
                                                              len(new_formatted_gpr),
                                                              new_formatted_gpr))

    return new_formatted_gpr


def expand_gpr(gpr):
    from sympy import expand
    formatted_gpr = gpr2expr(gpr)
    simplified_formatted_gpr = sympy.to_dnf(formatted_gpr)
    expanded_gpr = expr2gpr(simplified_formatted_gpr)

    return expanded_gpr


def expr2gpr(simplified_formatted_gpr):
    # reverse everything
    # multipliers_re = re.compile(r'\b(\d+\*)')
    # gid_re_rev = re.compile(r'g_(\d+)_(\d+)_g')
    # new_formatted_gpr = multipliers_re.sub('', simplified_formatted_gpr)
    # We deactivate this for sympy parsing of expressions
    # new_formatted_gpr = gid_re_rev.sub(r'\1.\2', new_formatted_gpr)
    # new_formatted_gpr = multiple_replace(new_formatted_gpr,EXPR2GPR_SUBS_DICT)
    new_formatted_gpr = multiple_replace(str(simplified_formatted_gpr),EXPR2GPR_SUBS_DICT)
    return new_formatted_gpr


def gpr2expr(gpr):
    # Recognize gene accession ids and rename them so that they don't look like numbers
    # gid_re = re.compile(r'(\d+)\.(\d+)')
    # formatted_gpr = gid_re.sub(r'g_\1_\2_g', gpr)
    formatted_gpr = str(gpr)
    formatted_gpr = multiple_replace(formatted_gpr,GPR2EXPR_SUBS_DICT)
    return formatted_gpr