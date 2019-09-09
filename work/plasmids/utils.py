# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Models vectors in ETFL models

.. moduleauthor:: ETFL team

Utility functions
"""


def read_seq(filename):
    with open(filename,'r') as fid:
        all_lines = fid.read_lines()

    seq = ''.join([x for line in all_lines for x in line.split()
                   if not x.isidigit()])

    return seq