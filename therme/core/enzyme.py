# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: ETFL team

ME-related Enzyme subclasses and methods definition


"""

from ..optim.variables import EnzymeVariable, ForwardEnzyme, BackwardEnzyme
from cobra import Species, Metabolite, DictList
from Bio.SeqUtils import molecular_weight
from .macromolecule import Macromolecule


class Enzyme(Macromolecule):
    def __init__(self, id=None, kcat=None, kcat_fwd=None, kcat_bwd=None,
                 kdeg=None, *args, **kwargs):
        Macromolecule.__init__(self, id = id, kdeg=kdeg, *args, **kwargs)

        if kcat is not None:
            self.kcat_fwd = kcat
            self.kcat_bwd = kcat
        else:
            if kcat_bwd is None and kcat_fwd is None:
                raise Exception('Either kcat must be provided, or both '
                                'kcat_fwb and kcat_bwd.')
            self.kcat_fwd = kcat_fwd
            self.kcat_bwd = kcat_bwd

        self.composition = None
        self.complexation = None

    def init_variable(self, queue=False):
        """
        Attach an EnzymeVariable object to the Species. Needs to have the object
        attached to a model

        :return:
        """
        self._internal_variable = self.model.add_variable(EnzymeVariable,
                                                        self,
                                                        queue=queue)

    @property
    def molecular_weight(self):
        # /!\ stoichiometric coefficient is negative
        return sum(-1*v*p.molecular_weight
                for p,v in self.complexation.metabolites.items())



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
    def from_metabolite(met, gene_id=None):

        new = Peptide(id=met.id,
                      name = met.name,
                      gene_id=gene_id)
        return new


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
