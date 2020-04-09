# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

ME-related Enzyme subclasses and methods definition


"""


from ..optim.variables import DNAVariable
from .macromolecule import Macromolecule

from Bio.SeqUtils import molecular_weight

class DNA(Macromolecule):
    def __init__(self, dna_len, gc_ratio, id='DNA', kdeg=0, *args, **kwargs):
        Macromolecule.__init__(self, id = id, kdeg=kdeg, *args, **kwargs)
        self.gc_ratio = gc_ratio
        self.len = dna_len


    def init_variable(self, queue=False):
        """
        Attach a DNAVariable object to the Species. Needs to have the object
        attached to a model

        :return:
        """
        self._internal_variable = self.model.add_variable(kind=DNAVariable,
                                    hook=self.model,
                                    id_=self.id,
                                    lb=0,
                                    ub=1)

    @property
    def molecular_weight(self):
        """
        Calculates the moleular weight of DNA based on the DNA GC-content and
        length

        :return:
        """

        g = self.gc_ratio

        # DNA mass (BioPython has g.mol^-1, while we are in mmol)
        ma = molecular_weight('A', seq_type='DNA') / 1000  # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
        mt = molecular_weight('T', seq_type='DNA') / 1000  # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
        mc = molecular_weight('C', seq_type='DNA') / 1000  # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
        mg = molecular_weight('G', seq_type='DNA') / 1000  # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1

        #         [     g.mmol(bp)^-1        * mmol(bp)/mmol(dna) ] ^ -1
        return ((1 - g) * (ma + mt) + g * (mc + mg)) * self.len


