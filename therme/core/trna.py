# -*- coding: utf-8 -*-
"""
.. module:: thermome
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: ETFL team

ME-related tRNA subclasses and methods definition


"""

from .macromolecule import Macromolecule
from ..optim.variables import tRNAVariable

class tRNA(Macromolecule):
    def __init__(self, aminoacid_id, charged, *args, **kwargs):
        if charged:
            prefix = 'charged'
        else:
            prefix = 'uncharged'

        aaid = prefix + '_tRNA_' + aminoacid_id
        Macromolecule.__init__(self, id = aaid, kdeg=0, *args, **kwargs)

        self._aminoacid_id = aminoacid_id

    def init_variable(self, queue=False):
        """
        Attach a tRNAVariable object to the Species. Needs to have the object
        attached to a model

        :return:
        """
        self._internal_variable = self.model.add_variable(tRNAVariable,
                                                      hook = self.model,
                                                      id_=self.id,
                                                      queue=queue)

    @property
    def aminoacid(self):
        return self.model.metabolites.get_by_id(self._aminoacid_id)