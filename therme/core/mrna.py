# -*- coding: utf-8 -*-
"""
.. module:: thermome
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

ME-related Enzyme subclasses and methods definition


"""

from ..optim.variables import mRNAVariable
from cobra import Species, Metabolite, DictList


class mRNA(Species):
    def __init__(self, id=None, kdeg=None, gene=None,*args, **kwargs):
        Species.__init__(self, id = id, *args, **kwargs)

        self.kdeg = kdeg
        self.gene = gene

    @property
    def peptide(self):
        return self.gene.peptide


    def init_variable(self, queue=False):
        """
        Attach an mRNAVariable object to the Species. Needs to have the object
        attached to a model

        :return:
        """
        self._mrna_variable = self.model.add_variable(mRNAVariable,
                                                        self,
                                                        queue=queue)

    @property
    def variable(self):
        """
        For convenience in the equations of the constraints

        :return:
        """
        try:
            return self._mrna_variable.variable
        except AttributeError:
            self.model.logger.warning('''{} has no model attached - variable attribute
             is not available'''.format(self.id))