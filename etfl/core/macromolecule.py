# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

ME-related macromolecule subclasses and methods definition


"""

from cobra import Species, Metabolite, DictList
from abc import ABC, abstractmethod, abstractproperty


class Macromolecule(Species, ABC):
    def __init__(self, id=None, kdeg=0, scaling_factor = None, *args, **kwargs):
        Species.__init__(self, id = id, *args, **kwargs)

        self.kdeg = kdeg
        self.degradation = None
        self._scaling_factor = scaling_factor

    @abstractmethod
    def init_variable(self, queue=False):
        """
        Attach an EnzymeVariable object to the Species. Needs to have the object
        attached to a model

        :return:
        """


    @property
    def concentration(self):
        """
        Concentration variable of the macromolecule in the cell.
        :return:
        """
        return self.scaling_factor * self.scaled_concentration


    @property
    def scaled_concentration(self):
        """
        Scaled concentration (scaling_factor*conc). If the scaling factor is the
        molecular weight, then this is similar to the mass fraction of the
        macromolecule in the cell, in g/gDW.
        :return:
        """
        return self.variable

    @property
    def X(self):
        """
        Value of the concentration after optimization.
        :return:
        """
        return self.scaling_factor * self.scaled_X


    @property
    def scaled_X(self):
        """
        Value of the scaled concentration (mass ratio) after optimization.
        :return:
        """
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
        # Calculate it only once
        if self._scaling_factor is None:
            self._scaling_factor = 1/self.molecular_weight
        return self._scaling_factor