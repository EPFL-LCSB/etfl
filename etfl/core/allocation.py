# -*- coding:utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

Core for the ME-part

"""


import numpy as np

from collections import defaultdict
from Bio.SeqUtils import molecular_weight

from tqdm import tqdm

from .genes import ExpressedGene
from .dna import DNA
from .rna  import mRNA
from .enzyme import Enzyme, Peptide
from .reactions import EnzymaticReaction, ProteinComplexation, \
    TranslationReaction, TranscriptionReaction, DegradationReaction, DNAFormation
from .expression import build_trna_charging, enzymes_to_gpr_no_stoichiometry, \
    make_stoich_from_aa_sequence, make_stoich_from_nt_sequence, \
    degrade_peptide, degrade_mrna, _extract_trna_from_reaction
from ..optim.constraints import SOS1Constraint, InterpolationConstraint, \
    mRNADegradation, EnzymeDegradation
from ..optim.variables import BinaryActivator, InterpolationVariable, DNAVariable, \
    GrowthRate, GenericVariable

from pytfa.optim.utils import chunk_sum, symbol_sum

MRNA_WEIGHT_CONS_ID = 'mRNA_weight_definition'
PROT_WEIGHT_CONS_ID = 'prot_weight_definition'
DNA_WEIGHT_CONS_ID  = 'DNA_weight_definition'
MRNA_WEIGHT_VAR_ID = 'mrna_ggdw'
PROT_WEIGHT_VAR_ID = 'prot_ggdw'
DNA_WEIGHT_VAR_ID  = 'dna_ggdw'
DNA_FORMATION_RXN_ID = 'DNA_formation'

def add_dummy_expression(model, aa_ratios, dummy_gene, dummy_peptide, dummy_protein,
                          peptide_length):
    # Making the reactions and adding them to the model so that they can
    # reference genes
    gtp = model.essentials['gtp']
    gdp = model.essentials['gdp']
    h2o = model.essentials['h2o']
    h = model.essentials['h']
    dummy_translation = TranslationReaction(id=model._get_translation_name(dummy_peptide),
                                            name='Dummy Translation',
                                            gene_id=dummy_gene.id,
                                            enzymes=model.ribosome.values(),
                                            scaled=True)
    dummy_complexation = ProteinComplexation(id='dummy_complexation',
                                             name='Dummy Complexation',
                                             target=dummy_protein,
                                             scaled=True)
    model.add_reactions([dummy_translation, dummy_complexation])
    model.translation_reactions += [dummy_translation]
    model.complexation_reactions += [dummy_complexation]
    # Use the input ratios to make the stoichiometry
    translation_mets = {}
    for k, v in aa_ratios.items():
        the_met_id = model.aa_dict[k]
        the_charged_trna, the_uncharged_trna, _ = model.trna_dict[the_met_id]
        translation_mets[the_charged_trna] = -1 * v * peptide_length
        translation_mets[the_uncharged_trna] = 1 * v * peptide_length
    translation_mets[model.metabolites.get_by_id(gtp)] = -2 * peptide_length
    translation_mets[model.metabolites.get_by_id(h2o)] = -2 * peptide_length
    translation_mets[model.metabolites.get_by_id(gdp)] = 2 * peptide_length
    translation_mets[model.metabolites.get_by_id(h)] = 2 * peptide_length
    # Do not forget to extract the tRNAs from the stoichiometry, since they
    # get diluted
    _extract_trna_from_reaction(translation_mets, dummy_translation)
    dummy_translation.add_metabolites(translation_mets, rescale=True)
    dummy_translation.add_peptide(dummy_peptide)
    dummy_complexation.add_peptides({dummy_peptide: -1})
    dummy_protein.init_variable()
    # Finally add the degradation flux
    prot_deg_stoich = dict()
    for k, v in aa_ratios.items():
        the_met_id = model.aa_dict[k]
        prot_deg_stoich[the_met_id] = v * peptide_length
    prot_deg_stoich[h2o] = -1 * peptide_length
    model._make_degradation_reaction(deg_stoich=prot_deg_stoich,
                                    macromolecule=dummy_protein,
                                    kind=EnzymeDegradation,
                                    scaled=True)
    # Now we can add the mass balance constraint
    model.add_mass_balance_constraint(dummy_complexation, dummy_protein)


def add_dummy_protein(model, dummy_peptide, enzyme_kdeg):
    # Create a dummy protein made of the dummy peptide
    dummy_protein = Enzyme(id='dummy_enzyme',
                           name='Dummy Enzyme',
                           kcat=0,
                           composition=[dummy_peptide.id],
                           kdeg=enzyme_kdeg)
    model.add_enzymes([dummy_protein])
    return dummy_protein


def add_dummy_peptide(model, aa_ratios, dummy_gene, peptide_length):
    # Create a dummy peptide
    dummy_peptide = Peptide(id='dummy_peptide',
                            name='Dummy peptide',
                            gene_id=dummy_gene.id)
    aa_weights = [v * molecular_weight(k, 'protein') for k, v in aa_ratios.items()]
    dummy_peptide.molecular_weight = peptide_length * sum(
        aa_weights) / 1000  # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
    dummy_peptide._model = model
    model.peptides += [dummy_peptide]
    return dummy_peptide


def add_dummy_mrna(model, dummy_gene, mrna_kdeg, mrna_length, nt_ratios):
    h2o = model.essentials['h2o']
    h = model.essentials['h']
    ppi = model.essentials['ppi']

    # Create a dummy mRNA
    dummy_mrna = mRNA(id='dummy_gene',
                      name='dummy mRNA',
                      kdeg=mrna_kdeg,
                      gene_id=dummy_gene.id)
    nt_weights = [v * molecular_weight(k, 'RNA') for k, v in nt_ratios.items()]
    dummy_mrna.molecular_weight = mrna_length * sum(
        nt_weights) / 1000  # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
    model.add_mrnas([dummy_mrna], add_degradation=False)
    dummy_transcription = TranscriptionReaction(id=model._get_transcription_name(dummy_mrna.id),
                                                name='Dummy Transcription',
                                                gene_id=dummy_gene.id,
                                                enzymes=model.rnap.values(),
                                                scaled=True)
    model.add_reactions([dummy_transcription])
    model.transcription_reactions += [dummy_transcription]
    # Use the input ratios to make the stoichiometry
    transcription_mets = {
        model.metabolites.get_by_id(model.rna_nucleotides[k]):
            -1 * v * mrna_length
        for k, v in nt_ratios.items()
    }
    transcription_mets[ppi] = mrna_length
    dummy_transcription.add_metabolites(transcription_mets, rescale=True)
    # Add the degradation
    mrna_deg_stoich = {
        model.metabolites.get_by_id(model.rna_nucleotides_mp[k]):
            -1 * v * mrna_length
        for k, v in nt_ratios.items()
    }
    mrna_deg_stoich[h2o] = -1 * mrna_length
    mrna_deg_stoich[h] = 1 * mrna_length
    model._make_degradation_reaction(deg_stoich=mrna_deg_stoich,
                                    macromolecule=dummy_mrna,
                                    kind=mRNADegradation,
                                    scaled=True)
    model.add_mass_balance_constraint(dummy_transcription, dummy_mrna)

    return dummy_mrna


def add_interpolation_variables(model):
    lambdas = []
    for e in range(model.n_mu_bins):
        lambda_i = model.add_variable(kind=BinaryActivator,
                                     hook=model,
                                     id_=str(e),
                                     lb=0,
                                     ub=1
                                     )
        lambdas += [lambda_i]
    sos_expr = symbol_sum(lambdas)

    model.add_constraint(kind=SOS1Constraint,
                        hook=model,
                        id_='interpolation_integer_SOS1',
                        expr=sos_expr,
                        lb=1,
                        ub=1)

    ga_vars = model.get_ordered_ga_vars()
    # mu_integer is the fraction coefficient of mu/mu_max:
    # mu_integer = delta_0*2^0 + delta_1*2^1 + ... + delta_n*2^n
    the_mu_integer = symbol_sum([(2 ** i) * ga_i
                                 for i, ga_i in enumerate(ga_vars)])

    # We want to equate the mu_integer with the currently active lambda index
    # 0*lambda_0 + 1*lambda_1 + ... + n*lambda_n = mu_integer
    ic_expr = symbol_sum([e * l for e, l in enumerate(lambdas)]) - the_mu_integer

    model.add_constraint(kind=InterpolationConstraint,
                        hook=model,
                        id_='growth_activators__EQ__interpolation_integers',
                        expr=ic_expr,
                        lb=0,
                        ub=0)

    model.regenerate_variables()
    model.regenerate_constraints()


def add_protein_mass_requirement(model, mu_values, p_rel):
    """
    Adds protein synthesis requirement

    input of type:

    ..code ::

        mu_values=[ 0.6,        1.0,        1.5,        2.0,        2.5     ]
        p_rel   = [ 0.675676,   0.604651,   0.540416,   0.530421,   0.520231]

        # mu_values in [h^-1]
        # p_rel in [g/gDw]

    :param mu_values:
    :param p_rel:
    :return:
    """

    activation_vars = model.get_variables_of_type(BinaryActivator)

    model_mus = [x[0] for x in model.mu_bins]
    p_hat = np.interp(x=model_mus,
                      xp=mu_values,
                      fp=p_rel)

    p_ref = symbol_sum([x * y for x, y in zip(p_hat, activation_vars)])

    # For legibility
    prot_ggdw = model.add_variable(kind=InterpolationVariable, hook=model,
                                  id_=PROT_WEIGHT_VAR_ID,
                                  lb=0,
                                  ub=1,  # can't have more prot than cell mass
                                  )

    epsilon = max(abs(np.diff(p_hat)))

    define_prot_weight_constraint(model, prot_ggdw)
    apply_prot_weight_constraint(model, p_ref, prot_ggdw, epsilon)

    model.interpolation_protein = p_hat
    model._interpolation_protein_tolerance = epsilon

    model.regenerate_variables()
    model.regenerate_constraints()


def apply_prot_weight_constraint(model, p_ref, prot_ggdw, epsilon):
    # E_ggdw = E_ref
    mass_coupling_expr = prot_ggdw - p_ref
    model.add_constraint(kind=InterpolationConstraint,
                         hook=model,
                         id_='prot_interpolation',
                         expr=mass_coupling_expr,
                         lb=-1 * epsilon,
                         ub=epsilon,
                         )


def define_prot_weight_constraint(model, prot_ggdw):
    # mmol.gDw^-1 / [scaling]
    enzyme_vars = model.enzymes.list_attr('concentration')
    # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
    enzyme_weights = model.enzymes.list_attr('molecular_weight')
    tot_prot = symbol_sum([x * y for x, y in zip(enzyme_weights, enzyme_vars)])
    # MW_1*[E1] + MW_2*[E2] + ... + MW_n*[En] = prot_ggdw
    mass_variable_def = tot_prot - prot_ggdw
    model.add_constraint(kind=InterpolationConstraint,
                        hook=model,
                        id_=PROT_WEIGHT_CONS_ID,
                        expr=mass_variable_def,
                        lb=0,
                        ub=0,
                        )


def add_rna_mass_requirement(model, mu_values, rna_rel):
    """
    Adds RNA synthesis requirement

    input of type:

    .. code::

        mu_values = [   0.6,        1.0,        1.5,        2.0,        2.5     ]
        rna_rel   = [   0.135135    0.151163    0.177829    0.205928    0.243931]

        # mu_values in [h^-1]
        # rna_rel in [g/gDw]

    :param mu_values:
    :param rna_rel:
    :return:
    """

    activation_vars = model.get_variables_of_type(BinaryActivator)

    model_mus = [x[0] for x in model.mu_bins]
    m_hat = np.interp(x=model_mus,
                      xp=mu_values,
                      fp=rna_rel)

    m_ref = symbol_sum([x * y for x, y in zip(m_hat, activation_vars)])

    # For legibility
    mrna_ggdw = model.add_variable(kind=InterpolationVariable,
                                  hook=model,
                                  id_=MRNA_WEIGHT_VAR_ID,
                                  lb=0,
                                  ub=1,  # can't have more rna than cell mass
                                  )

    epsilon = max(abs(np.diff(m_hat)))

    define_mrna_weight_constraint(model, mrna_ggdw)
    apply_mrna_weight_constraint(model, m_ref, mrna_ggdw, epsilon)

    model.interpolation_mrna = m_hat
    model._interpolation_mrna_tolerance = epsilon

    model.regenerate_variables()
    model.regenerate_constraints()


def apply_mrna_weight_constraint(model, m_ref, mrna_ggdw, epsilon):
    # mRNA_ggdw = mRNA_ref
    mass_coupling_expr = mrna_ggdw - m_ref
    model.add_constraint(kind=InterpolationConstraint,
                        hook=model,
                        id_='mRNA_interpolation',
                        expr=mass_coupling_expr,
                        lb=-1 * epsilon,
                        ub=epsilon,
                        )


def define_mrna_weight_constraint(model,mrna_ggdw):
    rna_vars = model.mrnas.list_attr('concentration')  # mmol.gDw^-1 / [scaling]
    rna_weights = model.mrnas.list_attr('molecular_weight')  # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
    tot_rna = symbol_sum([x * y for x, y in zip(rna_weights, rna_vars)])
    # MW_1*[rna1] + MW_2*[rna2] + ... + MW_n*[rna_n] = mRNA_ggdw
    mass_variable_def = tot_rna - mrna_ggdw
    model.add_constraint(kind=InterpolationConstraint,
                        hook=model,
                        id_=MRNA_WEIGHT_CONS_ID,
                        expr=mass_variable_def,
                        lb=0,
                        ub=0,
                        )


def add_dna_mass_requirement(model, mu_values, dna_rel, gc_ratio,
                             chromosome_len, dna_dict, ppi='ppi_c'):
    """
    Adds DNA synthesis requirement

    input of type:

    .. code::

        mu_values = [   0.6,        1.0,        1.5,        2.0,        2.5     ]
        dna_rel   = [   0.135135    0.151163    0.177829    0.205928    0.243931]

        # mu_values in [h^-1]
        # dna_rel in [g/gDw]

    :param mu_values:
    :param dna_rel:
    :return:
    """

    model.dna_nucleotides = dna_dict

    # Get mu interpolation
    activation_vars = model.get_variables_of_type(BinaryActivator)

    model_mus = [x[0] for x in model.mu_bins]
    m_hat = np.interp(x=model_mus,
                      xp=mu_values,
                      fp=dna_rel)

    m_ref = symbol_sum([x * y for x, y in zip(m_hat, activation_vars)])

    # Add DNA variable:

    dna = DNA(kdeg=0, dna_len=chromosome_len, gc_ratio=gc_ratio)
    # Assumption: kdeg for DNA is close to 0
    model.add_dna(dna)

    # Create dummy DNA reaction
    dna_formation = DNAFormation(id=DNA_FORMATION_RXN_ID, name='DNA Formation',
                                 dna=dna, mu_sigma= model._mu_range[-1],
                                 scaled=True)
    model.add_reactions([dna_formation])

    # In this formulation, we make 1 unit of the whole chromosome with NTPs
    g = gc_ratio
    mets = {v: -1 * chromosome_len * (g if k.lower() in 'gc' else 1 - g)
            for k, v in model.dna_nucleotides.items()}
    # Don't forget to release ppi (2 ppi per bp)
    mets[ppi] = 2 * chromosome_len

    dna_formation.add_metabolites(mets)

    # Add mass balance : 0 = v_syn - [mu]*[DNA]
    model.add_mass_balance_constraint(
        synthesis_flux=dna_formation,
        macromolecule=dna)

    # For legibility
    dna_ggdw = model.add_variable(kind=InterpolationVariable,
                                 hook=model,
                                 id_=DNA_WEIGHT_VAR_ID,
                                 lb=0,
                                 ub=1,  # can't have more dna than cell mass
                                 )

    define_dna_weight_constraint(model, dna, dna_ggdw, gc_ratio, chromosome_len)

    # DNA_ggdw = DNA_ref
    mass_coupling_expr = dna_ggdw - m_ref

    epsilon = max(abs(np.diff(m_hat)))

    model.add_constraint(kind=InterpolationConstraint,
                        hook=model,
                        id_='DNA_interpolation',
                        expr=mass_coupling_expr,
                        lb=-1 * epsilon,
                        ub=epsilon,
                        )

    model.interpolation_dna = m_hat
    model._interpolation_dna_tolerance = epsilon

    model.regenerate_variables()
    model.regenerate_constraints()


def define_dna_weight_constraint(model, dna, dna_ggdw, gc_content, chromosome_len):
    # DNA mass (BioPython has g.mol^-1, while we are in mmol)
    ma = molecular_weight('A', seq_type='DNA') / 1000  # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
    mt = molecular_weight('T', seq_type='DNA') / 1000  # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
    mc = molecular_weight('C', seq_type='DNA') / 1000  # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
    mg = molecular_weight('G', seq_type='DNA') / 1000  # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
    #              g.mmol(bp)^-1        * mmol(bp)/mmol(dna) * mmol(dna).gDW^-1
    tot_dna = ((1 - gc_content) * (ma + mt) + gc_content * (mc + mg)) * chromosome_len * dna.concentration
    # MW_avg*[DNA] = mRNA_ggdw
    # 1/scaling because the [X]s are scaled (eg mmol.ggDW^-1 -> back to mol.ggDW^1)
    mass_variable_def = tot_dna - dna_ggdw
    model.add_constraint(kind=InterpolationConstraint,
                         hook=model,
                         id_=DNA_WEIGHT_CONS_ID,
                         expr=mass_variable_def,
                         lb=0,
                         ub=0,
                         )