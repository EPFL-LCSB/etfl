# -*- coding: utf-8 -*-
"""
.. module:: thermome
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

ME-related Enzyme subclasses and methods definition


"""


from ..optim.variables import DNAVariable
from cobra import Species


class DNA(Species):
    def __init__(self, id='DNA', kdeg=0, *args, **kwargs):
        Species.__init__(self, id = id, *args, **kwargs)

        self.kdeg = kdeg


    def init_variable(self, queue=False):
        """
        Attach a DNAVariable object to the Species. Needs to have the object
        attached to a model

        :return:
        """
        self._dna_variable = self.model.add_variable(kind=DNAVariable,
                                    hook=self.model,
                                    id_=self.id,
                                    lb=0,
                                    ub=self.model.max_enzyme_concentration)

    @property
    def variable(self):
        """
        For convenience in the equations of the constraints

        :return:
        """
        try:
            return self._dna_variable.variable
        except AttributeError:
            self.throw_nomodel_error()


    def throw_nomodel_error(self):
        self.model.logger.warning('''{} has no model attached - variable attribute
         is not available'''.format(self.id))