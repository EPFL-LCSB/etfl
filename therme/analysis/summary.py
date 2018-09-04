# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Expression, Thermodynamics and Flux models

.. moduleauthor:: ETFL team

Summarizes quantities in models

"""

import pandas as pd

def get_amino_acid_consumption(model, solution = None,
                               trna_reaction_prefix='trna_ch_'):

    if solution is None:
        try:
            solution = model.solution
        except AttributeError as e:
            model.logger.error('If no solution object is provided, the model '
                               'must contain a solution attribute. You can get '
                               'on by running model.optimize()')
            raise(e)

    aminoacid_ids = list( model.aa_dict.values() )

    trna_charging_rxn_ids = [trna_reaction_prefix + x for x in aminoacid_ids]

    consumptions = [solution.x_dict[rid] for rid in trna_charging_rxn_ids]

    return pd.DataFrame(consumptions, index = aminoacid_ids, columns = ['tRNA_flux'])

def get_ntp_consumption(mode, solution = None):

    usage = dict()

    for ntp in model.rna_nucleotides:
         met = model.metabolites.get_by_id(ntp)
         usage[ntp] = sum(x.metabolites[met]*solution.flux[x.id]
                          for x in met.reactions
                          if x.id in model.transcription_reactions)

