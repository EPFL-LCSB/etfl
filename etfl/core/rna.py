# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

ME-related Enzyme subclasses and methods definition


"""
from cobra import Metabolite
from .macromolecule import Macromolecule
from ..optim.variables import mRNAVariable, tRNAVariable
from Bio.SeqUtils import molecular_weight
from .macromolecule import Macromolecule

from warnings import warn

class RNA(Macromolecule):
    def __init__(self, id=None, kdeg=None, gene_id=None, *args, **kwargs):
        Macromolecule.__init__(self, id = id, kdeg=kdeg, *args, **kwargs)

        self._gene_id = gene_id
        self._molecular_weight_override = 0

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
        self._internal_variable = \
            self.model.add_variable(mRNAVariable,
                                    self,
                                    scaling_factor=self.scaling_factor,
                                    lb = 0,
                                    ub=1,
                                    queue=queue)

    @property
    def molecular_weight(self):
        if not self._molecular_weight_override:
            return molecular_weight(self.rna, seq_type='RNA') / 1000 # g.mol^-1 ->
            # kg.mol^-1 (SI) =
            # g.mmol^-1
        else:
            return self._molecular_weight_override

    @molecular_weight.setter
    def molecular_weight(self, value):
        self._molecular_weight_override = value

class mRNA(RNA):
    # This class includes also rRNA and sRNA, because their mass balances are similar
    @property
    def peptide(self):
        return self.gene.peptide


class rRNA(Metabolite):
    def __init__(self, id=None, ribosomes=[], **kwargs):
        Metabolite.__init__(self, id=id, **kwargs)
        self._ribosomes = ribosomes
        
    @property
    def ribosomes(self):
        return self._ribosomes
    
    @ribosomes.setter
    def ribosomes(self,value):
        self._ribosomes = value
    
    @staticmethod
    def from_metabolite(met):

        new = rRNA(id=met.id,
                   name = met.name)
        new._model = met.model
        return new


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
        self._internal_variable = \
            self.model.add_variable(tRNAVariable,
                                    hook = self.model,
                                    id_=self.id,
                                    scaling_factor=self.scaling_factor,
                                    lb = 0,
                                    ub=1,
                                    queue=queue)

    @property
    def aminoacid(self):
        return self.model.metabolites.get_by_id(self._aminoacid_id)

    @property
    def molecular_weight(self):
        return 25 # an average weight based on https://bionumbers.hms.harvard.edu/files/Nucleic%20Acids_Sizes_and_Molecular_Weights_2pgs.pdf