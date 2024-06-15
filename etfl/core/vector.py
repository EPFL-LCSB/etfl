# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Models vectors in ETFL models

.. moduleauthor:: ETFL team

Model RNAP limitation in case of vector addition to a host
"""

from .rna import mRNA

from Bio.Seq import Seq

from argparse import ArgumentError

SEQ_TYPE_ERROR = "Vector type not recognized. Should be among " \
                 "plasmid, viral_dna, viral_rna"

def check_seq_type(sequence, vector_type):
    """
    Checks the sequence matches the type of the vector

    Plasmids are dsDNA,
    Viral vectors are DNA or RNA

    :param sequence:
    :param vector_type:
    :return:
    """
    l_vector_type = vector_type.lower()

    if isinstance(sequence, str):
        typed_sequence = make_seq(sequence, l_vector_type)
    elif isinstance(sequence, Seq):
        typed_sequence = sequence
    else:
        raise TypeError('The type of the sequence argument should be either '
                        'string or a Bio.Seq')

    # if l_vector_type == 'plasmid':
    #     assert(isinstance(sequence.alphabet, DNAAlphabet))
    # elif 'viral' in l_vector_type and 'rna' in l_vector_type:
    #     assert(isinstance(sequence.alphabet, RNAAlphabet))
    # elif 'viral' in l_vector_type and 'dna' in l_vector_type:
    #     assert(isinstance(sequence.alphabet, DNAAlphabet))
    # else:
    #     raise ArgumentError(SEQ_TYPE_ERROR)


    return typed_sequence


def make_seq(sequence, l_vector_type):
    if l_vector_type == 'plasmid':
        typed_sequence = Seq(sequence)#, alphabet=DNAAlphabet)
    elif 'viral' in l_vector_type and 'rna' in l_vector_type:
        typed_sequence = Seq(sequence)#, alphabet=RNAAlphabet)
    elif 'viral' in l_vector_type and 'dna' in l_vector_type:
        typed_sequence = Seq(sequence)#, alphabet=DNAAlphabet)
    else:
        raise ArgumentError(SEQ_TYPE_ERROR)
    return typed_sequence


def make_plasmid(gene_list):
    """
    gene_dict can be ordered to have an ordered sequence

    :param gene_list:
    :return:
    """
    sequence = ''.join([x.sequence * x.copy_number for x in gene_list])
    return sequence



class Vector:
    def __init__(self, id_, genes, reactions, gc_ratio, length,
                 proteins = None,
                 sequence=None,
                 mrna_dict=None,
                 coupling_dict=None,
                 rnap=None,
                 ribosome=None):
        """
        TODO: Explain in details

        :param sequence:
        :param genes:
        :param reactions:
        :param mrna_dict:
        :param coupling_dict:
        :param rnap:
        :param ribosome:
        """
        self.id = id_
        self.sequence = sequence
        self.genes = genes
        self.rnap = rnap
        self.ribosome = ribosome

        self.gc_ratio = gc_ratio
        self.len = length
        self.reactions = reactions
        self.mrna_dict = mrna_dict if mrna_dict is not None else dict()
        self.coupling_dict = coupling_dict if coupling_dict is not None else dict()
        self.proteins  = proteins if proteins is not None else list()
        
        self._peptides = []
        self._default_rnap = None
        self._default_rib = None
        
        self._calibration_tcpt = dict()
        self._calibration_tnsl = dict()

    def check_sequence(self):
        """
        To be defined in subclasses, check wehther RNA, DNA vector
        :return:
        """
        raise NotImplementedError
        
    def add_dna(self, dna):
        """
        Adds a DNA object to the model

        :param dna:
        :type dna: DNA
        :return:
        """

        try:
            dna._model = self.model
        except AttributeError:
            raise Exception('To define DNA species for the vector,'
                            'it should be associated to a model.')
        dna.init_variable()

        self.dna = dna

    @property
    def peptides(self):
        return self._peptides
    
    @peptides.setter
    def peptides(self, pep_list):
        if not hasattr(pep_list, '__iter__'):
            pep_list = [pep_list]
        self._peptides = pep_list
        
    @property
    def calibration_tcpt(self):
        return self._calibration_tcpt
    
    @calibration_tcpt.setter
    def calibration_tcpt(self, value):
        self._calibration_tcpt = value
        
    @property
    def calibration_tnsl(self):
        return self._calibration_tnsl
    
    @calibration_tnsl.setter
    def calibration_tnsl(self, value):
        self._calibration_tnsl = value
    
    @property
    def default_rnap(self):
        return self._default_rnap
    
    @default_rnap.setter
    def default_rnap(self, value):
        self._default_rnap = value
        
    @property
    def default_rib(self):
        return self._default_rib
    
    @default_rib.setter
    def default_rib(self, value):
        self._default_rib = value

    def build_default_mrna(self, kdeg, ids=None, force=False):
        """
        Build the mrna_dict attribute of the vector using a given value for
        the degradation rate.

        param kdeg: The mRNA degradation rate in h^-1
        param ids: list of genes for which to build mRNA objects. If None, all of them are used.
        param force: Allow overriding of existing mRNAs
        :return:
        """

        if ids is None:
            ids = [g.id for g in self.genes]

        for this_id in ids:
            self._add_mrna_to_dict(this_id, kdeg, force)

    def _add_mrna_to_dict(self, gene_id, kdeg, force):

        if gene_id in self.mrna_dict and not force:
            raise RuntimeError('There is already a mrna_dict entry for this gene'
                               'in the vector. Use the argument force = True to '
                               'override it.')

        mrna = mRNA(id=gene_id,
                    kdeg=kdeg,
                    gene_id=gene_id)

        self.mrna_dict[gene_id] = mrna



class Plasmid(Vector):
    def check_sequence(self, sequence):
        return check_seq_type(sequence, vector_type='plasmid')

    def integrate(self):
        pass

class ViralDNA(Vector):
    def check_sequence(self, sequence):
        return check_seq_type(sequence, vector_type='viral_dna')

    def integrate(self):
        pass

class ViralRNA(Vector):
    def check_sequence(self, sequence):
        return check_seq_type(sequence, vector_type='viral_rna')

    def integrate(self):
        pass
