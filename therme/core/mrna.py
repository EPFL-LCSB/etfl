# -*- coding: utf-8 -*-
"""
.. module:: thermome
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: ETFL team

ME-related Enzyme subclasses and methods definition


"""

from ..optim.variables import mRNAVariable
from Bio.SeqUtils import molecular_weight
from .macromolecule import Macromolecule



class mRNA(Macromolecule):
    def __init__(self, id=None, kdeg=None, gene_id=None, *args, **kwargs):
        Macromolecule.__init__(self, id = id, kdeg=kdeg, *args, **kwargs)

        self._gene_id = gene_id
        self._molecular_weight_override = 0

    @property
    def peptide(self):
        return self.gene.peptide

    @property
    def rna(self):
        return self.gene.rna

    @property
    def gene(self):
        return self.model.genes.get_by_id(self._gene_id)


    def init_variable(self, queue=False):
        """
        Attach an mRNAVariable object to the Species. Needs to have the object
        attached to a model

        :return:
        """
        self._internal_variable = self.model.add_variable(mRNAVariable,
                                                        self,
                                                        queue=queue)

    @property
    def molecular_weight(self):
        if not self._molecular_weight_override:
            return molecular_weight(self.rna) / 1000 # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
        else:
            return self._molecular_weight_override

    @molecular_weight.setter
    def molecular_weight(self, value):
        self._molecular_weight_override = value
