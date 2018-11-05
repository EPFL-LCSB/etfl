# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: ETFL team

ME-related macromolecule subclasses and methods definition


"""

from cobra import Species, Metabolite, DictList
from abc import ABC, abstractmethod, abstractproperty


class Macromolecule(Species, ABC):
    def __init__(self, id=None, kdeg=0, *args, **kwargs):
        Species.__init__(self, id = id, *args, **kwargs)

        self.kdeg = kdeg
        self.degradation = None

    @abstractmethod
    def init_variable(self, queue=False):
        """
        Attach an EnzymeVariable object to the Species. Needs to have the object
        attached to a model

        :return:
        """


    @property
    def concentration(self):
        return self.scaling_factor * self.scaled_concentration


    @property
    def scaled_concentration(self):
        return self.variable

    @property
    def X(self):
        return self.scaling_factor * self.scaled_X


    @property
    def scaled_X(self):
        return self.variable.primal

    @property
    def variable(self):
        """
        For convenience in the equations of the constraints

        :return:
        """
        try:
            return self._internal_variable.variable
        except AttributeError:
            self.throw_nomodel_error()

    def throw_nomodel_error(self):
        self.model.logger.warning('''{} has no model attached - variable attribute
         is not available'''.format(self.id))

    @abstractproperty
    def molecular_weight(self):
        """
        Necessary for scaling
        Use Biopython for this

        :return:
        """

    @property
    def scaling_factor(self):
        return 1/self.molecular_weight