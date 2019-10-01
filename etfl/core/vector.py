# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Models vectors in ETFL models

.. moduleauthor:: ETFL team

Model RNAP limitation in case of vector addition to a host
"""

from .memodel import MEModel
from .rna import mRNA
from .allocation import MRNA_WEIGHT_CONS_ID, PROT_WEIGHT_CONS_ID, \
    DNA_WEIGHT_CONS_ID, MRNA_WEIGHT_VAR_ID, PROT_WEIGHT_VAR_ID, \
    DNA_WEIGHT_VAR_ID, DNA_FORMATION_RXN_ID,\
    define_prot_weight_constraint, define_mrna_weight_constraint, \
    define_dna_weight_constraint
from ..optim.constraints import tRNAMassBalance, InterpolationConstraint
from ..optim.variables import InterpolationVariable

from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet, RNAAlphabet


SEQ_TYPE_ERROR = "Vector type ot recognized. Should be among " \
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

    if l_vector_type == 'plasmid':
        assert(isinstance(sequence.alphabet, DNAAlphabet))
    elif 'viral' in l_vector_type and 'rna' in l_vector_type:
        assert(isinstance(sequence.alphabet, RNAAlphabet))
    elif 'viral' in l_vector_type and 'dna' in l_vector_type:
        assert(isinstance(sequence.alphabet, DNAAlphabet))
    else:
        raise ArgumentError(SEQ_TYPE_ERROR)


    return typed_sequence


def make_seq(sequence, l_vector_type):
    if l_vector_type == 'plasmid':
        typed_sequence = Seq(sequence, alphabet=DNAAlphabet)
    elif 'viral' in l_vector_type and 'rna' in l_vector_type:
        typed_sequence = Seq(sequence, alphabet=RNAAlphabet)
    elif 'viral' in l_vector_type and 'dna' in l_vector_type:
        typed_sequence = Seq(sequence, alphabet=DNAAlphabet)
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

class TransModel(MEModel):

    def __init__(self, me_model, inplace = True):
        if not inplace:
            new = me_model.copy()
        else:
            new = me_model

        self._me_model = new
        self._has_transcription_changed = False
        self._has_translation_changed = False

    def __getattr__(self,attr):
        """
        Hack to subclass instantiated object:
        https://stackoverflow.com/questions/33463232/subclassing-in-python-of-instantiated-superclass

        :param attr:
        :return:
        """
        return getattr(self._me_model,attr)

    def add_vector(self, vector,
                   copy_number=1):

        self.add_reactions(vector.reactions)

        # Add the ribosome and RNAP of the vector if they exist
        if vector.rnap is not None:
            self.add_vector_RNAP(vector.rnap)
        if vector.ribosome is not None:
            self.add_vector_ribosome(vecor.ribosome)

        # Edit the gene copy number by the number of plasmids.
        for g in vector.genes:
            g.copy_number *= copy_number

        self.add_genes(vector.genes)
        self.express_genes(vector.genes)

        mrna_dict = vector.mrna_dict

        self.add_mrnas(mrna_dict.values())

        self.add_enzymatic_coupling(vector.coupling_dict)

        # recompute mRNA-dependent constraints
        self.recompute_trna_balances

        if self._has_transcription_changed:
            self.recompute_transcription()
        if self._has_translation_changed:
            self.recompute_translation()

        # This needs to account for new DNA
        self.recompute_allocation()

    def add_vector_RNAP(self, rnap):
        """
        Adds the vector's RNAP to the RNAP pool of the cell

        :param rnap:
        :return:
        """
        self.add_rnap(rnap, free_ratio=free_rnap_ratio)

        self._has_transcription_changed = True

        pass

    def add_vector_ribosome(self, ribosome):
        """
        Adds the vector's ribosome to the ribosome pool of the cell

        :param rnap:
        :return:
        """
        self.add_ribosome(ribosome, free_ratio=free_rib_ratio)

        self._has_translation_changed = True

        pass

    def recompute_transcription(self):
        """

        :return:
        """

        #1.1 Remove translation catalytic constraints

        #1.2 Apply new translation catalytic constraints

        #2.1 Remove former ribosome balance

        #2.2 Apply new ribosome balance

        pass

    def recompute_translation(self):
        """

        :return:
        """

        #1.1 Remove transcription catalytic constraints

        #1.2 Apply new transcription catalytic constraints

        #2.1 Remove former RNAP balance

        #2.2 Apply new RNAP balance
        pass

    def recompute_allocation(self):
        """

        :return:
        """

        #0. Does the model have allocation constraints ?

        #1. Remove previous allocation constraints

        #2. Apply new allocation constraints
        self.add_protein_mass_requirement(neidhardt_mu, neidhardt_prel)
        self.add_rna_mass_requirement(neidhardt_mu, neidhardt_rrel)
        self.add_dna_mass_requirement(mu_values=neidhardt_mu,
                                       dna_rel=neidhardt_drel,
                                       gc_ratio=gc_ratio,
                                       chromosome_len=chromosome_len,
                                       dna_dict=dna_nucleotides)


class Vector:
    def __init__(self, sequence, genes, reactions,
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
        self.sequence = self.check_sequence(sequence)
        self.genes = genes
        self.rnap = rnap
        self.ribosome = ribosome

        self.reactions = reactions
        self.mrna_dict = mrna_dict if mrna_dict is not None else dict()
        self.coupling_dict = coupling_dict if coupling_dict is not None else dict()

    def check_sequence(self):
        """
        To be defined in subclasses, check wehther RNA, DNA vector
        :return:
        """
        raise NotImplementedError

    def integrate(self):
        """
        DNA vectors are added to chromosome length
        RNA Vectors are added to RNA mass

        :return:
        """
        raise NotImplementedError

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
            self._add_mrna_to_dict(this_id, kdeg)

    def _add_mrna_to_dict(self, gene_id, kdeg):

        if this_id in self.mrna_dict and not force:
            raise RuntimeError('There is already a mrna_dict entry for this gene'
                               'in the vector. Use the argument force = True to '
                               'override it.')

        mrna = mRNA(id=gene_id,
                    kdeg=kdeg,
                    gene_id=gene_id)

        self.mrna_dict[gene_id] = mrna



class Plasmid(Vector):
    def check_sequence(self, sequence):
        return check_seq_type(sequence, type='plasmid')

    def integrate(self):
        pass

class ViralDNA(Vector):
    def check_sequence(self, sequence):
        return check_seq_type(sequence, type='viral_dna')

    def integrate(self):
        pass

class ViralRNA(Vector):
    def check_sequence(self, sequence):
        return check_seq_type(sequence, type='viral_rna')

    def integrate(self):
        pass
