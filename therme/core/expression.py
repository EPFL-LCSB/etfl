# -*- coding: utf-8 -*-
"""
.. module:: thermome
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

ME-related Reaction subclasses and methods definition


"""
from cobra import Reaction, Metabolite

from collections import defaultdict

def build_trna_charging(model, aa_dict,
                        atp='atp_c',
                        amp='amp_c',
                        ppi='ppi_c',
                        h2o='h2o_c',
                        h='h_c'):
    trna_dict = dict()
    for letter, aa_id in aa_dict.items():
        aa = model.metabolites.get_by_id(aa_id)

        rxn_id = 'trna_ch_{}'.format(aa.id)

        charged_trna = Metabolite(name = 'Charged tRNA-{}'.format(aa.name),
                                  id = 'trna_charged_{}'.format(aa.id),
                                  compartment = aa.compartment)
        uncharged_trna = Metabolite(name = 'Uncharged tRNA-{}'.format(aa.name),
                                  id = 'trna_uncharged_{}'.format(aa.id),
                                  compartment = aa.compartment)
        charging_rxn = Reaction(name='tRNA Charging of {}'.format(aa.name),
                                id= rxn_id)

        model.add_reactions([charging_rxn])

        # Trick in case two amino acids are linked to the same reaction. Example:
        # Cysteine and selenocysteine
        the_rxn = model.reactions.get_by_id(rxn_id)

        the_rxn.add_metabolites({
            aa:-1,
            uncharged_trna:-1,
            atp:-1,
            h2o:-2,
            charged_trna:1,
            amp:1,
            ppi:1,
            h:2,
        })
        trna_dict[aa] = (charged_trna,uncharged_trna)
    return trna_dict

def make_stoich_from_aa_sequence(sequence, model, aa_dict, trna_dict,
                                 gtp, gdp, h2o, h):
    stoich = defaultdict(int)

    for letter in sequence:
        met_id = aa_dict[letter]
        met = model.metabolites.get_by_id(met_id)
        charged_trna, uncharged_trna = trna_dict[met]
        # stoich[met]-=1
        stoich[charged_trna] -= 1
        stoich[uncharged_trna] += 1
    stoich[gtp] = -2 * len(sequence)
    stoich[h2o] = -2 * len(sequence)
    stoich[gdp] = 2 * len(sequence)
    stoich[h] = 2 * len(sequence)
    return stoich

def make_stoich_from_nt_sequence(sequence, model, nt_dict):
    stoich = defaultdict(int)
    for letter in sequence:
        met_id = nt_dict[letter]
        met = model.metabolites.get_by_id(met_id)
        stoich[met]-=1
    return stoich