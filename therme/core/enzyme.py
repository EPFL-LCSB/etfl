# -*- coding: utf-8 -*-
"""
.. module:: thermome
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

ME-related Enzyme subclasses and methods definition


"""

from ..optim.variables import EnzymeVariable, ForwardEnzyme, BackwardEnzyme
from ..optim.constraints import TotalEnzyme
from cobra import Species, Metabolite, DictList
from Bio.SeqUtils import molecular_weight


class Enzyme(Species):
    def __init__(self, id=None, kcat=None, kcat_fwd=None, kcat_bwd=None,
                 kdeg=None, *args, **kwargs):
        Species.__init__(self, id = id, *args, **kwargs)

        if kcat is not None:
            self.kcat_fwd = kcat
            self.kcat_bwd = kcat
        else:
            if kcat_bwd is None and kcat_fwd is None:
                raise Exception('Either kcat must be provided, or both '
                                'kcat_fwb and kcat_bwd.')
            self.kcat_fwd = kcat_fwd
            self.kcat_bwd = kcat_bwd

        self.kdeg = kdeg
        self.composition = None
        self.complexation = None


    def init_variable(self, queue=False):
        """
        Attach an EnzymeVariable object to the Species. Needs to have the object
        attached to a model

        :return:
        """
        self._enzyme_variable = self.model.add_variable(EnzymeVariable,
                                                        self,
                                                        queue=False)

    @property
    def variable(self):
        """
        For convenience in the equations of the constraints

        :return:
        """
        try:
            return self._enzyme_variable.variable
        except AttributeError:
            self.throw_nomodel_error()

    @property
    def molecular_weight(self):
        # /!\ stoichiometric coefficient is negative
        return sum(-1*v*p.molecular_weight
                for p,v in self.complexation.metabolites.items())


    def throw_nomodel_error(self):
        self.model.logger.warning('''{} has no model attached - variable attribute
         is not available'''.format(self.id))




class Peptide(Metabolite):
    """
    Subclass to describe peptides resulting from gene translation
    """

    def __init__(self, id=None, gene_id=None, **kwargs):
        Metabolite.__init__(self, id=id, **kwargs)
        self._gene_id = gene_id
        self._molecular_weight_override = 0


    @property
    def gene(self):
        try:
            return self.model.genes.get_by_id(self._gene_id)
        except KeyError:
            self.model.logger.warning('Peptide {} tried to reference {}, '
                                      'not in model'.format(self.id, self._gene_id))
            return None

    @property
    def peptide(self):
        return self.gene.peptide

    @property
    def molecular_weight(self):
        if not self._molecular_weight_override:
            return molecular_weight(self.peptide) / 1000 # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
        else:
            return self._molecular_weight_override

    @molecular_weight.setter
    def molecular_weight(self, value):
        self._molecular_weight_override = value

    @staticmethod
    def from_metabolite(met, copy = False):
        if copy:
            met = met.copy()
        met.__class__ = Peptide
        return met


class Ribosome(Enzyme):
    def __init__(self, id=None, kribo=None, kdeg=None, *args, **kwargs):
        Enzyme.__init__(self, id = id, kdeg=kdeg, kcat = kribo, *args, **kwargs)

    @property
    def kribo(self):
        return self.kcat_fwd


class RNAPolymerase(Enzyme):
    def __init__(self, id=None, ktrans=None, kdeg=None, *args, **kwargs):
        Enzyme.__init__(self, id = id, kdeg=kdeg, kcat = ktrans, *args, **kwargs)

    @property
    def ktrans(self):
        return self.kcat_fwd
