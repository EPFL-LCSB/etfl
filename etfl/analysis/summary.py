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

    solution = check_solution(model, solution)

    aminoacid_ids = list( model.aa_dict.values() )

    trna_charging_rxn_ids = [trna_reaction_prefix + x for x in aminoacid_ids]

    consumptions = [solution.x_dict[rid] for rid in trna_charging_rxn_ids]

    return pd.DataFrame(consumptions, index = aminoacid_ids, columns = ['tRNA_flux'])

def get_ntp_consumption(model, solution = None):

    usage = dict()

    solution = check_solution(model, solution)

    ntp_ids = model.rna_nucleotides.values()
    for ntp in ntp_ids:
         met = model.metabolites.get_by_id(ntp)
         usage[ntp] = sum(x.metabolites[met]*solution.fluxes[x.id]
                          for x in met.reactions
                          if x.id in model.transcription_reactions)

         return pd.DataFrame(usage, index=ntp_ids, columns=['ntp_usage'])


def check_solution(model, solution):
    if solution is None:
        try:
            solution = model.solution
        except AttributeError as e:
            model.logger.error('If no solution object is provided, the model '
                               'must contain a solution attribute. You can get '
                               'on by running model.optimize()')
            raise (e)
    return solution


def print_standard_sol(model, solution = None, flux_dict = None):

    solution = check_solution(model, solution)

    if flux_dict is None:
        flux_dict = {'Growth': model.growth_reaction.id}
        if 'EX_glc__D_e' in model.reactions:
            flux_dict['Glucose uptake'] = 'EX_glc__D_e'

    try: width = max([max([len(x) for x in the_dict]) for the_dict in [flux_dict,
                                                                 model.ribosome,
                                                                 model.rnap ]])
    except ValueError: #Happens when dict empty:
        width = 10

    print('{0: <{width}}: {1}'.format('Objective', solution.objective_value,
                                      width=width+3))
    # print()
    # print('Fluxes')
    _print_dict_items_fluxes(solution, flux_dict, width)
    # print()
    # print('Ribosomes produced')
    _print_dict_items_vars(solution, model.ribosome, width)
    # print()
    # print('RNAP produced')
    _print_dict_items_vars(solution, model.rnap, width)
    try:
        # print()
        print(' - {0: >{width}}: {1}'.format('DNA produced',
                                          solution.raw[model.dna.variable.name],
                                          width=width))

    except AttributeError:
        pass


def _print_dict_items_vars(solution, the_dict, width):
    if the_dict:
        for k, v in the_dict.items():
            x = solution.raw[v.variable.name] * v.scaling_factor
            print(' - {0: >{width}}: {1}'.format(k, x, width=width))


def _print_dict_items_fluxes(solution, the_dict, width):
    if the_dict:
        for k, v in the_dict.items():
            x = solution.fluxes[v]
            print(' - {0: >{width}}: {1}'.format(k, x, width=width))