# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

ME-related Reaction subclasses and methods definition


"""
from warnings import warn
from cobra import Gene
from Bio.Seq import Seq
#from Bio.Alphabet import DNAAlphabet, RNAAlphabet, ProteinAlphabet
from ..optim.constraints import SynthesisConstraint, MinimalAllocation, MinimalCoupling


def make_sequence(sequence):#, seq_type):
    """
    seq_type must be an instance of DNAAphabet(), RNAAlphabet, or ProteinAlphabet
    :param sequence:
    :param seq_type:
    :return:
    """
    if isinstance(sequence, str):
        typed_sequence = Seq(sequence)#, seq_type())
    elif isinstance(sequence, Seq):
        #assert (isinstance(sequence.alphabet, seq_type))
        typed_sequence = sequence
    else:
        raise TypeError('The type of the sequence argument should be either '
                        'string or a Bio.Seq')
    
    return typed_sequence

class ExpressedGene(Gene):
    '''This calss represents the genes that can be transcribed'''
    def __init__(self, id, name, sequence, copy_number=1,
                 transcribed_by=None, min_tcpt_activity=0,
                 *args, **kwargs):
        Gene.__init__(self, id, name, *args, **kwargs)
        self.sequence = make_sequence(sequence)#, DNAAlphabet)
        self._rna = ''

        self._copy_number = int(copy_number)
        self._transcribed_by = transcribed_by
        self._min_tcpt_activity = min_tcpt_activity

    @property
    def copy_number(self):
        return self._copy_number

    @copy_number.setter
    def copy_number(self, value):
        # TODO: Make this a setter that rewrites the adequate constraints
        if value != self._copy_number:
            if self.model is None:
                # Easy
                self._copy_number = int(value)
            else:
                # We need to make the model change the RNAP allocation
                self._copy_number = int(value)
                self.model.edit_gene_copy_number(self.id)
        else:
            # Nothing to do here :)
            pass

    @property
    def transcribed_by(self):
        return self._transcribed_by

    @transcribed_by.setter
    def transcribed_by(self,value):
        if value != self._transcribed_by:
            if self.model is None:
                # Easy
                self._transcribed_by = value
            else:
                mod_id = self.id + '_transcription'
                cstr = self.model.get_constraints_of_type(SynthesisConstraint)
                cstr_exist = True # an indicator to show if the gene is transcribed
                try:
                    cstr.get_by_id(mod_id)
                except KeyError:
                    cstr_exist = False
                if cstr_exist:
                    # trancription constraint already exists
                    # TODO: We need to make the model change the corresponding constraints
                    raise NotImplementedError()
                else:
                    self._transcribed_by = value
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
        self._rna=make_sequence(value)#,RNAAlphabet)
        
    @property
    def min_tcpt_activity(self):
        return self._min_tcpt_activity

    @min_tcpt_activity.setter
    def min_tcpt_activity(self,value):
        if value != self._min_tcpt_activity:
            if self.model is None:
                # Easy
                self._min_tcpt_activity = value
            else:
                mod_id = self.id
                cstr = self.model.get_constraints_of_type(MinimalAllocation)
                cstr_exist = True # an indicator to show if the gene is transcribed
                try:
                    cstr.get_by_id(mod_id)
                except KeyError:
                    cstr_exist = False
                if cstr_exist:
                    # minimal allocation constraint already exists
                    warn('The minimal allocation constraint might not be changed.'
                         'To change it, use the change_tcpt_basal_activity function.')
                    self._min_tcpt_activity = value
                else:
                    self._min_tcpt_activity = value
        else:
            # Nothing to do here :)
            pass
    
    
class CodingGene(ExpressedGene):
    '''This calss represents the genes that can be translated into protein'''
    def __init__(self, id, name, sequence, min_tnsl_activity=0,
                 translated_by=None,
                 *args, **kwargs):
        ExpressedGene.__init__(self, id, name, sequence, *args, **kwargs)

        self._peptide = ''
        self._translated_by = translated_by
        self._min_tnsl_activity = min_tnsl_activity
    
    @property
    def translated_by(self):
        return self._translated_by

    @translated_by.setter
    def translated_by(self,value):
        if value != self._translated_by:
            if self.model is None:
                # Easy
                self._translated_by = value
            else:
                # We need to make the model change the corresponding cstr
                mod_id = self.id + '_translation'
                cstr = self.model.get_constraints_of_type(SynthesisConstraint)
                cstr_exist = True # an indicator to show if the gene is translated
                try:
                    cstr.get_by_id(mod_id)
                except KeyError:
                    cstr_exist = False
                if cstr_exist:
                    # translation constraint already exists
                    # TODO: We need to make the model change the corresponding constraints
                    raise NotImplementedError()
                else:
                    self._translated_by = value
        else:
            # Nothing to do here :)
            pass
        
    @property
    def peptide(self):
        if not self._peptide:
            # Translation table 11 is for bacteria
            the_pep = self.rna.translate(to_stop = False, table = 11)
            the_pep = str(the_pep).replace('*','')
            self._peptide = make_sequence(the_pep)#, ProteinAlphabet)

        return self._peptide

    @peptide.setter
    def peptide(self,value):
        self._peptide=make_sequence(value)#,ProteinAlphabet)
        
    @property
    def min_tnsl_activity(self):
        return self._min_tnsl_activity

    @min_tnsl_activity.setter
    def min_tnsl_activity(self,value):
        if value != self._min_tnsl_activity:
            if self.model is None:
                # Easy
                self._min_tnsl_activity = value
            else:
                mod_id = self.id
                cstr = self.model.get_constraints_of_type(MinimalCoupling)
                cstr_exist = True # an indicator to show if the gene is transcribed
                try:
                    cstr.get_by_id(mod_id)
                except KeyError:
                    cstr_exist = False
                if cstr_exist:
                    # minimal coupling constraint already exists
                    warn('The minimal coupling constraint might not be changed.'
                         'To change it, use the change_tnsl_basal_activity function.')
                    self._min_tnsl_activity = value
                else:
                    self._min_tnsl_activity = value
        else:
            # Nothing to do here :)
            pass
        
    @staticmethod
    def from_gene(gene, sequence):
        """
        This method clones a cobra.Gene object into an CodingGene,
        and attaches a sequence to it

        :param gene: the gene to reproduce
        :type gene: cobra.Gene
        :param sequence: a string-like dna sequence
        :return: an CodingGene  object
        """
        new =  CodingGene(   id = gene.id,
                                name= gene.name,
                                sequence=sequence,
                                functional = gene.functional)

        return new
        