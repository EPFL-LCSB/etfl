# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Adding vector to an ETFL model

.. moduleauthor:: ETFL team

Model RNAP limitation in case of vector addition to a host
"""

from .memodel import MEModel
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet, RNAAlphabet


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

    if l_vector_type == 'plasmid':
        typed_sequence = Seq(sequence, alphabet = DNAAlphabet)
    elif 'viral' in l_vector_type and 'rna' in l_vector_type:
        typed_sequence = Seq(sequence, alphabet = RNAAlphabet)
    elif 'viral' in l_vector_type and 'dna' in l_vector_type:
        typed_sequence = Seq(sequence, alphabet = DNAAlphabet)
    else:
        raise ArgumentError("Vector type ot recognized. Should be among "
                            "plasmid, viral_dna, viral_rna")


    return typed_sequence

def make_plasmid(gene_dict):
    """
    gene_dict can be ordered to have an ordered sequence

    :param gene_dict:
    :return:
    """
    sequence = ''.join([x.sequence * x.copy_number for x in gene_dict.values()])
    return sequence

class TransModel(MEModel):

    def __init__(self, me_model, inplace = True):
        if not inplace:
            new = me_model.copy()
        else:
            new = me_model

        self._me_model = new

    def __getattr__(self,attr):
        """
        Hack to subclass instantiated object:
        https://stackoverflow.com/questions/33463232/subclassing-in-python-of-instantiated-superclass

        :param attr:
        :return:
        """
        return getattr(self._me_model,attr)

    def add_vector(self, gene_dict, sequence, vector_type=None, copy_number=1):

        if vector_type is not None:
            sequence = check_seq_type(sequence, vector_type)

        pass

    def add_vector_RNAP(self, rnap):
        """
        Adds the vector's RNAP to the RNAP pool of the cell

        :param rnap:
        :return:
        """

        pass

    def add_vector_ribosome(self, ribosome):
        """
        Adds the vector's ribosome to the ribosome pool of the cell

        :param rnap:
        :return:
        """

        pass


