# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: ETFL team

ME-related Enzyme subclasses and methods definition


"""


from ..optim.variables import DNAVariable
from .macromolecule import Macromolecule


class DNA(Macromolecule):
    def __init__(self, id='DNA', kdeg=0, *args, **kwargs):
        Macromolecule.__init__(self, id = id, kdeg=kdeg, *args, **kwargs)


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
                                    ub=self.model.max_enzyme_concentration)

