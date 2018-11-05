# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: ETFL team

ME-related Reaction subclasses and methods definition


"""

from cobra import Gene
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet, ProteinAlphabet

class ExpressedGene(Gene):
    """
    Subclass to describe reactions that are catalyzed by an enzyme.
    """
    def __init__(self, id, name, sequence, *args, **kwargs):
        Gene.__init__(self, id, name, *args, **kwargs)
        self.sequence = Seq(sequence, DNAAlphabet())
        self._rna = ''
        self._peptide = ''

    @property
    def rna(self):
        if not self._rna:
            self._rna = self.sequence.transcribe()

        return self._rna

    @property
    def peptide(self):
        if not self._peptide:
            # Translation table 11 is for bacteria
            the_pep = self.rna.translate(to_stop = False, table = 11)
            the_pep = str(the_pep).replace('*','')
            self._peptide = Seq(the_pep, ProteinAlphabet())

        return self._peptide

    @staticmethod
    def from_gene(gene, sequence):
        """
        This method clones a cobra.Gene object into an ExpressedGene,
        and attaches a sequence to it

        :param gene: the gene to reproduce
        :type gene: cobra.Gene
        :param sequence: a string-like dna sequence
        :return: an Expressed Gene  object
        """
        new =  ExpressedGene(   id = gene.id,
                                name= gene.name,
                                sequence=sequence,
                                functional = gene.functional)

        return new