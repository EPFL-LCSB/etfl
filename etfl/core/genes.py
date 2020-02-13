# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

ME-related Reaction subclasses and methods definition


"""

from cobra import Gene
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet, RNAAlphabet, ProteinAlphabet
from etfl.optim.constraints import SynthesisConstraint


def make_sequence(sequence, seq_type):
    """
    seq_type must be an instance of DNAAphabet(), RNAAlphabet, or ProteinAlphabet
    :param sequence:
    :param seq_type:
    :return:
    """
    if isinstance(sequence, str):
        typed_sequence = Seq(sequence, seq_type())
    elif isinstance(sequence, Seq):
        assert (isinstance(sequence.alphabet, seq_type))
        typed_sequence = sequence
    else:
        raise TypeError('The type of the sequence argument should be either '
                        'string or a Bio.Seq')
    
    return typed_sequence

class ExpressedGene(Gene):
    """
    Subclass to describe reactions that are catalyzed by an enzyme.
    """
    def __init__(self, id, name, sequence, copy_number=1,
                 transcribed_by=None,
                 translated_by=None,
                 *args, **kwargs):
        Gene.__init__(self, id, name, *args, **kwargs)
        self.sequence = make_sequence(sequence, DNAAlphabet)
        self._rna = ''
        self._peptide = ''

        self._copy_number = copy_number
        self._transcribed_by = transcribed_by
        self._translated_by = translated_by

    @property
    def copy_number(self):
        return self._copy_number

    @copy_number.setter
    def copy_number(self, value):
        # TODO: Make this a setter that rewrites the adequate constraints
        if value != self._copy_number:
            if self.model is None:
                # Easy
                self._copy_number = value
            else:
                # We need to make the model change the RNAP allocation
                self._copy_number = value
                self.model.edit_gene_copy_number(self.id)
        else:
            # Nothing to do here :)
            pass

    @property
    def transcribed_by(self):
        return self._transcribed_by

    @transcribed_by.setter
    def transcribed_by(self,value):
        # TODO: Make this a setter that rewrites the adequate constraints
        if value != self._transcribed_by:
            if self.model is None:
                # Easy
                self._transcribed_by = value
            else:
                # We need to make the model change the corresponding cstr
                mod_id = self.id + '_transcription'
                Cstr = self.model.get_constraints_of_type(SynthesisConstraint)
                cstr_exist = 1
                try:
                    Cstr.get_by_id(mod_id)
                except KeyError:
                    cstr_exist = 0
                if cstr_exist:
                    # trancription constraint already exists
                    raise NotImplementedError()
                else:
                    self._transcribed_by = value
        else:
            # Nothing to do here :)
            pass

    @property
    def translated_by(self):
        return self._translated_by

    @translated_by.setter
    def translated_by(self,value):
        # TODO: Make this a setter that rewrites the adequate constraints
        if value != self._translated_by:
            if self.model is None:
                # Easy
                self._translated_by = value
            else:
                # We need to make the model change the corresponding cstr
                mod_id = self.id + '_translation'
                Cstr = self.model.get_constraints_of_type(SynthesisConstraint)
                cstr_exist = 1
                try:
                    Cstr.get_by_id(mod_id)
                except KeyError:
                    cstr_exist = 0
                if cstr_exist:
                    # trancription constraint already exists
                    raise NotImplementedError()
                else:
                    self._translated_by = value
        else:
            # Nothing to do here :)
            pass


    @property
    def rna(self):
        if not self._rna:
            self._rna = self.sequence.transcribe()

        return self._rna

    @rna.setter
    def rna(self,value):
        self._rna=make_sequence(value,RNAAlphabet)

    @property
    def peptide(self):
        if not self._peptide:
            # Translation table 11 is for bacteria
            the_pep = self.rna.translate(to_stop = False, table = 11)
            the_pep = str(the_pep).replace('*','')
            self._peptide = make_sequence(the_pep, ProteinAlphabet)

        return self._peptide

    @peptide.setter
    def peptide(self,value):
        self._peptide=make_sequence(value,ProteinAlphabet)

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