# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: ETFL team

Analysis utilities

"""

from collections import defaultdict

def enzymes_to_peptides_conc(model, enzyme_conc):
    """

    :param enzyme_conc: dict or pandas.Series, with the key/index being enzyme
                        variable names, and the value their concentration
    :return:
    """

    peptide_conc = defaultdict(float)

    for enzyme_var_name, conc in enzyme_conc.items():
        the_var = model._var_dict[enzyme_var_name]
        the_enzyme = model.enzymes.get_by_id(the_var.id)
        complex_dict = the_enzyme.complexation.metabolites

        for pep, stoich in complex_dict.items():
            peptide_conc[pep.id] += abs(stoich) * conc

    return peptide_conc