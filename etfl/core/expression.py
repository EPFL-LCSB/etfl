# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

ME-related Reaction subclasses and methods definition


"""
from cobra import Reaction, Metabolite

from ..core.rna import tRNA
from .rna import tRNA
from ..utils.parsing import parse_gpr

from warnings import warn
from collections import defaultdict

def build_trna_charging(model, aa_dict,
                        atp='atp_c',
                        amp='amp_c',
                        ppi='ppi_c',
                        h2o='h2o_c',
                        h='h_c'):
    """
    Build th tRNA charging reactions, based on the amino acid dictionary

    :param model: An ETFL Model
    :type model: :class:`etfl.core.memodel.MEModel`
    :param aa_dict: A dictionary of aminoacid letter to amicoacid met id

        **Example :**

        .. code-block:: python

            aa_dict = {
                        'A':'ala__L_c',
                        'R':'arg__L_c',
                        ...
                    }
    
    :param atp: metabolite ID of the cytosolic ATP
    :param amp: metabolite ID of the cytosolic AMP
    :param ppi: metabolite ID of the cytosolic diphosphate
    :param h2o: metabolite ID of the cytosolic water
    :param h: metabolite ID of the cytosolic hydrogen ions
    :return: A dictionary of tRNAs, keys are the aminoacido letters,
            values the charging reactions
    """
    trna_dict = dict()
    for letter, aa_id in aa_dict.items():
        aa = model.metabolites.get_by_id(aa_id)

        rxn_id = get_trna_charging_id(aa_id)

        # charged_trna = Metabolite(name = 'Charged tRNA-{}'.format(aa.name),
        #                           id = 'trna_charged_{}'.format(aa.id),
        #                           compartment = aa.compartment)
        # uncharged_trna = Metabolite(name = 'Uncharged tRNA-{}'.format(aa.name),
        #                           id = 'trna_uncharged_{}'.format(aa.id),
        #                           compartment = aa.compartment)
        charging_rxn = Reaction(name='tRNA Charging of {}'.format(aa.name),
                                id= rxn_id)

        charged_trna = tRNA(aminoacid_id=aa.id,
                            charged = True,
                            name = aa.name)

        uncharged_trna = tRNA(aminoacid_id=aa.id,
                            charged = False,
                            name = aa.name)

        model.add_reactions([charging_rxn])

        # Trick in case two amino acids are linked to the same reaction. Example:
        # Cysteine and selenocysteine
        the_rxn = model.reactions.get_by_id(rxn_id)
        mets = {
                aa:-1,
                # uncharged_trna:-1,
                atp:-1,
                h2o:-2,
                # charged_trna:1,
                amp:1,
                ppi:1,
                h:2,
                }

        the_rxn.add_metabolites(mets)

        trna_dict[aa_id] = (charged_trna,uncharged_trna, charging_rxn)
    return trna_dict


def get_trna_charging_id(aa_id):
    rxn_id = 'trna_ch_{}'.format(aa_id)
    return rxn_id


def make_stoich_from_aa_sequence(sequence, aa_dict, trna_dict,
                                 gtp, gdp, pi, h2o, h):
    """
    Makes the stoichiometry of the peptide synthesis reaction based on the
    amino acid sequence

    :param sequence: sequence of aminoacids (letter form)
    :type sequence: :class:`Bio.Seq` or :class:`str`
    :param aa_dict: A dictionary of aminoacid letter to amicoacid met id

        **Example :**

        .. code-block:: python

            aa_dict = {
                        'A':'ala__L_c',
                        'R':'arg__L_c',
                        ...
                    }
    
    :param trna_dict: the dict returned by :func:`etfl.core.expression.build_trna_charging`
    :param gtp: metabolite ID for GTP
    :param gdp: metabolite ID for GDP
    :param pi: metabolite ID for phosphate
    :param h2o: metabolite ID for water
    :param h: metabolite ID for H+
    :return:
    """
    stoich = defaultdict(int)

    for letter in sequence:
        met_id = aa_dict[letter]
        charged_trna, uncharged_trna, _ = trna_dict[met_id]
        # stoich[met]-=1
        stoich[charged_trna] -= 1
        stoich[uncharged_trna] += 1
    stoich[gtp] = -2 * len(sequence)
    stoich[h2o] = -2 * len(sequence)
    stoich[gdp] = 2 * len(sequence)
    stoich[h] = 2 * len(sequence)
    stoich[pi] = 2 * len(sequence)
    return stoich

def make_stoich_from_nt_sequence(sequence, nt_dict, ppi):
    """
    Makes the stoichiometry of the RNA synthesis reaction based on the
    nucleotides sequence

    :param sequence: sequence of RNA nucleotides
    :type sequence: :class:`Bio.Seq` or :class:`str`
    :param nt_dict: A dictionary of RNA nucleotide triphosphate
                            letter to nucleotideTP met id
        **Example :**

        .. code-block:: python

            rna_nucleotides = {
                        'A':'atp_c',
                        'U':'utp_c',
                        ...
                    }
    
    :param ppi: metabolite ID for diphosphate
    :return:
    """
    stoich = defaultdict(int)
    for letter in sequence:
        met_id = nt_dict[letter]
        stoich[met_id]-=1
    stoich[ppi] = len(sequence)
    return stoich

def degrade_peptide(peptide, aa_dict, h2o):
    """
    Degrades a peptide in amino acids, based on its sequence

    :param peptide: The peptide
    :type peptide: :class:`etfl.core.enzyme.Peptide`
    :param aa_dict: A dictionary of aminoacid letter to amicoacid met id

        ** Example : **

        .. code-block:: python

            aa_dict = {
                        'A':'ala__L_c',
                        'R':'arg__L_c',
                        ...
                    }

    :param h2o: metabolite ID for water
    :return:
    """
    sequence = peptide.peptide

    stoich = defaultdict(int)

    for letter in sequence:
        met_id = aa_dict[letter]
        stoich[met_id]+=1
    stoich[h2o] = -1 * len(sequence)

    return stoich

def degrade_mrna(mrna, nt_dict, h2o, h):
    """
    Degrades a mRNA in nucleotides monophosphate, based on its sequence

    :param mrna: The peptide
    :type mrna: :class:`etfl.core.rna.mRNA`
    :param nt_dict: A dictionary of RNA nucleotide monophosphate letter to nucleotideMP met id

        **Example :**

        .. code-block:: python

            rna_nucleotides_mp = {
                        'A':'amp_c',
                        'U':'ump_c',
                        ...
                    }

    :param h2o: metabolite ID for water
    :param h: metabolite ID for H+
    :return:
    """
    sequence = mrna.rna

    stoich = defaultdict(int)

    for letter in sequence:
        met_id = nt_dict[letter]
        stoich[met_id]+=1
    stoich[h2o] = -1 * len(sequence)
    stoich[h] = 1 * len(sequence)

    return stoich


def is_me_compatible(reaction):
    """
    Check if a Cobra reaction has sufficient information to add expression coupling

    :param reaction:
    :type reaction: :class:`cobra.core.Reaction`
    :return:
    """
    # Test if the GPR is a proper one:
    this_gpr = reaction.gene_reaction_rule
    is_proper_gpr = bool(this_gpr) and this_gpr != '[]' and not 's0001' in this_gpr

    # sym_gpr = parse_gpr(this_gpr)

    ret = True

    if not is_proper_gpr:
        # Then we cannot constrain
        warn('Improper GPR for {}'.format(reaction.id))
        ret = False

    # # Check that all the genes participating in this gpr have a translation
    # # reaction:
    # is_translated = {x: '{}_translation'.format(x.name) \
    #                     in translation_reactions
    #                  for x in sym_gpr.free_symbols}
    # if not all(is_translated.values()):
    #     warn(
    #         'Not all peptides in the GPR of {} are translated: {}'.format(
    #             reaction.id, is_translated))
    #     ret = False

    return ret

def enzymes_to_gpr(rxn):
    """
    Builds a Gene to Protein to Reaction association rules from the enzymes of
    an enzymatic reaction

    :param rxn:
    :return:
    """
    return ' or '.join([' ( ' + ' and '.join(['{}*{}'.format(v,k)
                                for v,k in isozyme.composition.items()])  + ' ) '
                                for isozyme in rxn.enzymes] )

def enzymes_to_gpr_no_stoichiometry(rxn):
    """
    Builds a Gene to Protein to Reaction association rules from the enzymes of
    an enzymatic reaction

    :param rxn:
    :return:
    """
    return ' or '.join([' ( ' + ' and '.join([v
                                for v in isozyme.composition])  + ' ) '
                                for isozyme in rxn.enzymes] )


def _extract_trna_from_reaction(aa_stoichiometry, rxn):
    """
    Read a stoichiometry dictionary, and replaces free aminoacids with tRNAs

    :param aa_stoichiometry: the stoichiometry dict to edit
    :type aa_stoichiometry: (dict) {:class:`cobra.core.Metabolite`: Number}
    :param rxn: the reaction whose stoichiometry is inspected
    :type rxn: :class:`cobra.core.Reaction`
    :return:
    """
    # Extract the tRNAs, since they will be used for a different mass balance
    # in self.add_trna_mass_balances
    for met, stoich in list(aa_stoichiometry.items()):
        if isinstance(met, tRNA):
            rxn.trna_stoich[met.id] = aa_stoichiometry.pop(met)