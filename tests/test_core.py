# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

Utilities to create a small model from a 1 reaction model in FBA

"""

from etfl.tests.small_model import create_etfl_model


def test_etfl():
    create_etfl_model(has_thermo=False,has_neidhardt=True)