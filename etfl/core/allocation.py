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
from cobra import Metabolite, Reaction

from tqdm import tqdm

from .carbohydrate import Carbohydrate
from .ion import Ion
from .lipid import Lipid
from .dna import DNA
from .rna  import mRNA
from .enzyme import Enzyme, Peptide
from .reactions import EnzymaticReaction, ProteinComplexation, \
    TranslationReaction, TranscriptionReaction, DegradationReaction, DNAFormation
from .expression import build_trna_charging, enzymes_to_gpr_no_stoichiometry, \
    make_stoich_from_aa_sequence, make_stoich_from_nt_sequence, \
    degrade_peptide, degrade_mrna, _extract_trna_from_reaction
from ..optim.constraints import SOS1Constraint, InterpolationConstraint, \
    mRNADegradation, EnzymeDegradation, ConstantAllocation
from ..optim.variables import BinaryActivator, InterpolationVariable, DNAVariable, \
    GrowthRate, GenericVariable, EnzymeVariable, mRNAVariable

from pytfa.optim.utils import chunk_sum, symbol_sum

PROT_CONSTANT_CONS_ID = 'prot_fix'
RNA_CONSTANT_CONS_ID = 'RNA_fix'
DNA_CONSTANT_CONS_ID = 'DNA_fix'
ENZYME_CONSTANT_CONS_ID = 'enzyme_fix'
TRNA_CONSTANT_CONS_ID = 'trna_fix'
MRNA_WEIGHT_CONS_ID = 'mRNA_weight_definition'
PROT_WEIGHT_CONS_ID = 'prot_weight_definition'
DNA_WEIGHT_CONS_ID  = 'DNA_weight_definition'
MRNA_WEIGHT_VAR_ID = 'mrna_ggdw'
PROT_WEIGHT_VAR_ID = 'prot_ggdw'
DNA_WEIGHT_VAR_ID  = 'dna_ggdw'
DNA_FORMATION_RXN_ID = 'DNA_formation'
LIPID_FORMATION_RXN_ID = 'Lipid_formation'
LIPID_WEIGHT_VAR_ID  = 'lipid_ggdw'
LIPID_WEIGHT_CONS_ID = 'lipid_weight_definition'
ION_FORMATION_RXN_ID = 'ion_formation'
ION_WEIGHT_VAR_ID  = 'ion_ggdw'
ION_WEIGHT_CONS_ID = 'ion_weight_definition'
CARBOHYDRATE_FORMATION_RXN_ID = 'Carbohydrate_formation'
CARBOHYDRATE_WEIGHT_VAR_ID  = 'carbohydrate_ggdw'
CARBOHYDRATE_WEIGHT_CONS_ID = 'carbohydrate_weight_definition'
VECTOR_CONSTANT_CONS_ID = 'vector_fix'
VECTOR_FORMATION_RXN_ID = 'vector_formation'

def fix_prot_ratio(model, prot_ratio):
    '''
    To keep consistency between FBA and ETFL biomass compositions, we divide biomass
    into two parts: BS1 and BS2. BS1 includes variable parts of biomass (i.e. RNA
    and protein), while BS2 includes the other components that are not modeled
    explicitly.
    inputs:
        model: ME-model
        mass_ratios: a dict of mass_ratios for biomass composition in the GEM
            It must have ratios for 'RNA' and 'protein'. If 'total mass' is provided,
            it is used to scale ratios. Otherwise, it's assumed to be 1 gr.
    outputs:
        return a model with an additional constraint on sum of RNA and protein share
    '''
        
                
    # mmol.gDw^-1 / [scaling]
    enzyme_vars = model.enzymes.list_attr('concentration')
    # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
    enzyme_weights = model.enzymes.list_attr('molecular_weight')
    
    # In ribsosome molecular weight the rRNAs are included but they are not protein, so should be removed
    ribosomal_rrna = find_rrna_weight_in_rib(model)
    
    expr = symbol_sum([x * y for x, y in zip(enzyme_weights, enzyme_vars)]) - ribosomal_rrna
    
    
    model.add_constraint(kind = ConstantAllocation, 
                                 hook = model, 
                                 expr = expr,
                                 id_ = PROT_CONSTANT_CONS_ID,
                                 lb = prot_ratio,
                                 ub = prot_ratio)
    
def fix_RNA_ratio(model, rna_ratio):
    '''
    To keep consistency between FBA and ETFL biomass compositions, we divide biomass
    into two parts: BS1 and BS2. BS1 includes variable parts of biomass (i.e. RNA
    and protein), while BS2 includes the other components that are not modeled
    explicitly.
    inputs:
        model: ME-model
        mass_ratios: a dict of mass_ratios for biomass composition in the GEM
            It must have ratios for 'RNA' and 'protein'. If 'total mass' is provided,
            it is used to scale ratios. Otherwise, it's assumed to be 1 gr.
    outputs:
        return a model with an additional constraint on sum of RNA and protein share
    '''
        
                
    rna_vars = model.mrnas.list_attr('concentration')  # mmol.gDw^-1 / [scaling]
    rna_weights = model.mrnas.list_attr('molecular_weight')  # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
    expr = symbol_sum([x * y for x, y in zip(rna_weights, rna_vars)])
    
    # part of rRNA is used to produce ribosomes; this should be considered
    ribosomal_rrna = find_rrna_weight_in_rib(model)
        
    expr = expr + ribosomal_rrna
    
    model.add_constraint(kind = ConstantAllocation, 
                                 hook = model, 
                                 expr = expr,
                                 id_ = RNA_CONSTANT_CONS_ID,
                                 lb = rna_ratio,
                                 ub = rna_ratio)
    
def fix_DNA_ratio(model, dna_ratio, gc_ratio, chromosome_len, tol = 0.05):
    '''
    A function similar  to fix_RNA_ratio. Used only in the case of adding vector
    and when variable biomass composition is not available. It adds a DNA species
    to the model that with a constant concentration, but this can be used for
    RNAP allocation constraints (to be compatible with those constraints).
    tol: a tolerance ration for the deviation of DNA from its mass ratio
    '''
        
    
    dna = DNA(kdeg=0, dna_len=chromosome_len, gc_ratio=gc_ratio)
    # Assumption: kdeg for DNA is close to 0
    model.add_dna(dna)          
    
    model.add_constraint(kind = ConstantAllocation, 
                                 hook = model, 
                                 expr = dna.variable,
                                 id_ = DNA_CONSTANT_CONS_ID,
                                 lb = dna_ratio - tol * dna_ratio,
                                 ub = dna_ratio + tol * dna_ratio)

def constrain_enzymes(model, enz_ratio, prot_ratio=None):
    # a function to add a constraint on total amount of enzymes based on their
    # fraction from total amount of proteins (should be before adding dummy)
                
    enz_vars = model.get_variables_of_type(EnzymeVariable)
    
    # we should first exclude dummy, ribosomes and rnaps
    exclusion = ['dummy_enzyme', #'rib', 'rib_mit', 'rnap', 'rnap_mit'
                 ]
    exclusion = ['EZ_{}'.format(x) for x in exclusion]
    enz_vars = [x for x in enz_vars if x.name not in exclusion]
    
    # In ribsosome molecular weight the rRNAs are included but they are not protein, so should be removed
    ribosomal_rrna = find_rrna_weight_in_rib(model)
    
    expr = symbol_sum([x for x in enz_vars]) - ribosomal_rrna
    try:
        # model with variable biomass
        model.add_constraint(kind = ConstantAllocation, 
                                 hook = model, 
                                 expr = expr - enz_ratio * model.variables.IV_prot_ggdw,
                                 id_ = ENZYME_CONSTANT_CONS_ID,
                                 ub = 0)
    except AttributeError:
        model.add_constraint(kind = ConstantAllocation, 
                                 hook = model, 
                                 expr = expr,
                                 id_ = ENZYME_CONSTANT_CONS_ID,
                                 lb = 0, # cannot be negative
                                 ub = enz_ratio * prot_ratio)
        
def constrain_trna(model, trna_ratio, rna_ratio=None):
    # Since solely based on mass balances we do not have a non-zero concentraio
    # for tRNA, an empirical allocation constraints should be added
                
    trna_var = model.get_variables_of_type(mRNAVariable).get_by_id('tRNA_gene')
    
    
    try:
        # model with variable biomass
        model.add_constraint(kind = ConstantAllocation, 
                                 hook = model, 
                                 expr = trna_var - trna_ratio * model.variables.IV_mrna_ggdw,
                                 id_ = TRNA_CONSTANT_CONS_ID,
                                 lb = 0,
                                 ub = 0)
    except AttributeError:
        model.add_constraint(kind = ConstantAllocation, 
                                 hook = model, 
                                 expr = trna_var,
                                 id_ = TRNA_CONSTANT_CONS_ID,
                                 lb = trna_ratio * rna_ratio,
                                 ub = trna_ratio * rna_ratio)

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
    model.add_enzymes([dummy_protein], prep = False)
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


def add_dummy_mrna(model, dummy_gene, mrna_kdeg, mrna_length, nt_ratios, name=''):
    h2o = model.essentials['h2o']
    h = model.essentials['h']
    ppi = model.essentials['ppi']

    # Create a dummy mRNA
    dummy_mrna = mRNA(id=dummy_gene.id,
                      name='dummy {}'.format(name),
                      kdeg=mrna_kdeg,
                      gene_id=dummy_gene.id)
    nt_weights = [v * molecular_weight(k, 'RNA') for k, v in nt_ratios.items()]
    dummy_mrna.molecular_weight = mrna_length * sum(
        nt_weights) / 1000  # g.mol^-1 -> kg.mol^-1 (SI) = g.mmol^-1
    model.add_mrnas([dummy_mrna], add_degradation=False)
    dummy_transcription = TranscriptionReaction(id=model._get_transcription_name(dummy_mrna.id),
                                                name='Dummy {} Transcription'.format(name),
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
    
     # In ribsosome molecular weight the rRNAs are included but they are not protein, so should be removed
    ribosomal_rrna = find_rrna_weight_in_rib(model)
        
    tot_prot = symbol_sum([x * y for x, y in zip(enzyme_weights, enzyme_vars)]) - ribosomal_rrna
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
    # part of rRNA is used to produce ribosomes; this should be considered
    ribosomal_rrna = find_rrna_weight_in_rib(model)
        
    tot_rna = symbol_sum([x * y for x, y in zip(rna_weights, rna_vars)]) + ribosomal_rrna
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

    mets = get_dna_synthesis_mets(model, chromosome_len, gc_ratio, ppi)

    dna_formation.add_metabolites(mets)

    # Add mass balance : 0 = v_syn - [mu]*[DNA]
    model.add_mass_balance_constraint(
        synthesis_flux=dna_formation,
        macromolecule=dna)
    
    epsilon = max(abs(np.diff(m_hat)))

    # For legibility
    dna_ggdw = model.add_variable(kind=InterpolationVariable,
                                 hook=model,
                                 id_=DNA_WEIGHT_VAR_ID,
                                 lb=0,
                                 ub=1,  # can't have more dna than cell mass
                                 )

    define_dna_weight_constraint(model, dna, dna_ggdw, gc_ratio, chromosome_len)
    apply_dna_weight_constraint(model, m_ref, dna_ggdw, epsilon)

    model.interpolation_dna = m_hat
    model._interpolation_dna_tolerance = epsilon

    model.regenerate_variables()
    model.regenerate_constraints()
    
def get_dna_synthesis_mets(model, chromosome_len, gc_ratio, ppi):
    # In this formulation, we make 1 unit of the whole chromosome with NTPs
    g = gc_ratio
    mets = {v: -1 * chromosome_len * (g if k.lower() in 'gc' else 1 - g)
            for k, v in model.dna_nucleotides.items()}
    # Don't forget to release ppi (2 ppi per bp)
    mets[ppi] = 2 * chromosome_len
    return mets
    
def apply_dna_weight_constraint(model, m_ref, dna_ggdw, epsilon):
    # DNA_ggdw = DNA_ref
    mass_coupling_expr = dna_ggdw - m_ref

    model.add_constraint(kind=InterpolationConstraint,
                        hook=model,
                        id_='DNA_interpolation',
                        expr=mass_coupling_expr,
                        lb=-1 * epsilon,
                        ub=epsilon,
                        )

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
    
def add_lipid_mass_requirement(model, lipid_mets, mass_ratio, mu_values,
                               lipid_rel, lipid_rxn = None):
    '''
    In general, we have two main situations:
        1) the lipid paripates in biomass formation as lumped metabolite.
        2) the lipid components partipate in biomass formation individually.
    In the first case, we should remove lipid metabolite from the model and replace
    it with a mcromolecule with a new mass balnce constraint.
    In the second case, after removing lipid metabolites from biomass rxn,
    we should define a new reaction to lump lipid metabolites. Then, it becomes
    similar to the first case.

    Parameters
    ----------
    model : MeModel
        ETFL model with variable biomass composition.
    lipid_mets : list
        A list of lipid metabolite id(s)
    mass_ratios : dict
        Keys are strings for biomass components and values are their ration in
        FBA model. The ratios should be consistent with the current stoichiometric
        coefficients. 
    mu_values : list or DataFrame
        Values of growth rates for which experimental data is available 
    lipid_rel : ist or DataFrame
        Different ratios of lipid for different growth rates
    lipid_rxn : string
        the rxn id for lipid psedoreaction. If None, there is no such reaction.

    Returns
    -------
    None.

    '''
    biomass_rxn = model.growth_reaction
    if isinstance(lipid_rxn, str): # the lumped lipid metabolite case
        lipid_id = lipid_mets[0]
        lipid = model.metabolites.get_by_id(lipid_id)
        # assumption: stoichiometric coefficient for lipid in biomass rxn is 1
        # to remove it from the model
        model.remove_metabolites(lipid)
        # to create a dict of lipid composition
        lipid_formation = model.reactions.get_by_id(lipid_rxn)
        # lipid_formation.id = LIPID_FORMATION_RXN_ID
        composition_dict = lipid_formation.metabolites
    else: # the individual lipid metabolites case
        mets = [model.metabolites.get_by_id(x) for x in lipid_mets]
        composition_dict = {met:biomass_rxn.get_coefficient(met) for met \
                      in mets}
        # to remove them from biomass rxn
        biomass_rxn.subtract_metabolites(composition_dict)
        # Define a dummy reation to lump components
        lipid_formation = Reaction(id = LIPID_FORMATION_RXN_ID,
                                   name = 'Lipid formation',
                                  lower_bound = - model.bigM,
                                  upper_bound = model.bigM)
        model.add_reactions([lipid_formation])
        lipid_formation.add_metabolites(composition_dict)
        
    lipid = Lipid(kdeg = 0, composition = composition_dict,
                  mass_ratio = mass_ratio)
    model.add_lipid(lipid)
    
    activation_vars = model.get_variables_of_type(BinaryActivator)
    model_mus = [x[0] for x in model.mu_bins]
    l_hat = np.interp(x=model_mus,
                      xp=mu_values,
                      fp=lipid_rel)

    l_ref = symbol_sum([x * y for x, y in zip(l_hat, activation_vars)])
    epsilon = max(abs(np.diff(l_hat)))
    
    
    model.add_mass_balance_constraint(synthesis_flux=lipid_formation,
                                      macromolecule=lipid)
    
    apply_lipid_weight_constraint(model, l_ref, lipid, epsilon)

    model.interpolation_lipid = l_hat
    model._interpolation_lipid_tolerance = epsilon

    model.regenerate_variables()
    model.regenerate_constraints()
    
def apply_lipid_weight_constraint(model, l_ref, lipid, epsilon):
    # mRNA_ggdw = mRNA_ref
    tot_lipid = lipid.variable
    mass_coupling_expr = tot_lipid - l_ref
    model.add_constraint(kind=InterpolationConstraint,
                        hook=model,
                        id_='lipid_interpolation',
                        expr=mass_coupling_expr,
                        lb=-1 * epsilon,
                        ub=epsilon,
                        )
    
def add_carbohydrate_mass_requirement(model, carbohydrate_mets, mass_ratio, mu_values,
                               carbohydrate_rel, carbohydrate_rxn = None):
    '''
    In general, we have two main situations:
        1) the carbohydrate paripates in biomass formation as lumped metabolite.
        2) the carbohydrate components partipate in biomass formation individually.
    In the first case, we should remove carbohydrate metabolite from the model and replace
    it with a mcromolecule with a new mass balnce constraint.
    In the second case, after removing carbohydrate metabolites from biomass rxn,
    we should define a new reaction to lump carbohydrate metabolites. Then, it becomes
    similar to the first case.

    Parameters
    ----------
    model : MeModel
        ETFL model with variable biomass composition.
    carbohydrate_mets : list
        A list of carbohydrate metabolite id(s)
    mass_ratios : dict
        Keys are strings for biomass components and values are their ration in
        FBA model. The ratios should be consistent with the current stoichiometric
        coefficients. 
    mu_values : list or DataFrame
        Values of growth rates for which experimental data is available 
    carbohydrate_rel : ist or DataFrame
        Different ratios of carbohydrate for different growth rates
    carbohydrate_rxn : string
        the rxn id for carbohydrate psedoreaction. If None, there is no such reaction.

    Returns
    -------
    None.

    '''
    biomass_rxn = model.growth_reaction
    if isinstance(carbohydrate_rxn, str): # the lumped carbohydrate metabolite case
        carbohydrate_id = carbohydrate_mets[0]
        carbohydrate = model.metabolites.get_by_id(carbohydrate_id)
        # assumption: stoichiometric coefficient for carbohydrate in biomass rxn is 1
        # to remove it from the model
        model.remove_metabolites(carbohydrate)
        # to create a dict of carbohydrate composition
        carbohydrate_formation = model.reactions.get_by_id(carbohydrate_rxn)
        # carbohydrate_formation.id = carbohydrate_FORMATION_RXN_ID
        composition_dict = carbohydrate_formation.metabolites
    else: # the individual carbohydrate metabolites case
        mets = [model.metabolites.get_by_id(x) for x in carbohydrate_mets]
        composition_dict = {met:biomass_rxn.get_coefficient(met) for met \
                      in mets}
        # to remove them from biomass rxn
        biomass_rxn.subtract_metabolites(composition_dict)
        # Define a dummy reation to lump components
        carbohydrate_formation = Reaction(id = CARBOHYDRATE_FORMATION_RXN_ID,
                                   name = 'Carbohydrate formation',
                                  lower_bound = - model.bigM,
                                  upper_bound = model.bigM)
        model.add_reactions([carbohydrate_formation])
        carbohydrate_formation.add_metabolites(composition_dict)
        
    carbohydrate = Carbohydrate(kdeg = 0, composition = composition_dict,
                  mass_ratio = mass_ratio)
    model.add_carbohydrate(carbohydrate)
    
    activation_vars = model.get_variables_of_type(BinaryActivator)
    model_mus = [x[0] for x in model.mu_bins]
    c_hat = np.interp(x=model_mus,
                      xp=mu_values,
                      fp=carbohydrate_rel)

    c_ref = symbol_sum([x * y for x, y in zip(c_hat, activation_vars)])
    epsilon = max(abs(np.diff(c_hat)))
    
    model.add_mass_balance_constraint(synthesis_flux=carbohydrate_formation,
                                      macromolecule=carbohydrate)
    
    apply_carbohydrate_weight_constraint(model, c_ref, carbohydrate, epsilon)

    model.interpolation_carbohydrate = c_hat
    model._interpolation_carbohydrate_tolerance = epsilon

    model.regenerate_variables()
    model.regenerate_constraints()
    
def apply_carbohydrate_weight_constraint(model, c_ref, carbohydrate, epsilon):
    # mRNA_ggdw = mRNA_ref
    tot_carbohydrate = carbohydrate.variable
    mass_coupling_expr = tot_carbohydrate - c_ref
    model.add_constraint(kind=InterpolationConstraint,
                        hook=model,
                        id_='carbohydrate_interpolation',
                        expr=mass_coupling_expr,
                        lb=-1 * epsilon,
                        ub=epsilon,
                        )
    
def add_ion_mass_requirement(model, ion_mets, mass_ratio, mu_values,
                               ion_rel, ion_rxn = None):
    '''
    In general, we have two main situations:
        1) the ion paripates in biomass formation as lumped metabolite.
        2) the ion components partipate in biomass formation individually.
    In the first case, we should remove ion metabolite from the model and replace
    it with a mcromolecule with a new mass balnce constraint.
    In the second case, after removing ion metabolites from biomass rxn,
    we should define a new reaction to lump ion metabolites. Then, it becomes
    similar to the first case.

    Parameters
    ----------
    model : MeModel
        ETFL model with variable biomass composition.
    ion_mets : list
        A list of ion metabolite id(s)
    mass_ratios : dict
        Keys are strings for biomass components and values are their ration in
        FBA model. The ratios should be consistent with the current stoichiometric
        coefficients. 
    mu_values : list or DataFrame
        Values of growth rates for which experimental data is available 
    ion_rel : ist or DataFrame
        Different ratios of ion for different growth rates
    ion_rxn : string
        the rxn id for ion psedoreaction. If None, there is no such reaction.

    Returns
    -------
    None.

    '''
    biomass_rxn = model.growth_reaction
    if isinstance(ion_rxn, str): # the lumped ion metabolite case
        ion_id = ion_mets[0]
        ion = model.metabolites.get_by_id(ion_id)
        # assumption: stoichiometric coefficient for ion in biomass rxn is 1
        # to remove it from the model
        model.remove_metabolites(ion)
        # to create a dict of ion composition
        ion_formation = model.reactions.get_by_id(ion_rxn)
        # ion_formation.id = ION_FORMATION_RXN_ID
        composition_dict = ion_formation.metabolites
    else: # the individual ion metabolites case
        mets = [model.metabolites.get_by_id(x) for x in ion_mets]
        composition_dict = {met:biomass_rxn.get_coefficient(met) for met \
                      in mets}
        # to remove them from biomass rxn
        biomass_rxn.subtract_metabolites(composition_dict)
        # Define a dummy reation to lump components
        ion_formation = Reaction(id = ION_FORMATION_RXN_ID,
                                   name = 'ion formation',
                                  lower_bound = - model.bigM,
                                  upper_bound = model.bigM)
        model.add_reactions([ion_formation])
        ion_formation.add_metabolites(composition_dict)
        
    ion = Ion(kdeg = 0, composition = composition_dict,
                  mass_ratio = mass_ratio)
    model.add_ion(ion)
    
    activation_vars = model.get_variables_of_type(BinaryActivator)
    model_mus = [x[0] for x in model.mu_bins]
    i_hat = np.interp(x=model_mus,
                      xp=mu_values,
                      fp=ion_rel)

    i_ref = symbol_sum([x * y for x, y in zip(i_hat, activation_vars)])
    epsilon = max(abs(np.diff(i_hat)))
    
    model.add_mass_balance_constraint(synthesis_flux=ion_formation,
                                      macromolecule=ion)
    
    apply_ion_weight_constraint(model, i_ref, ion, epsilon)

    model.interpolation_ion = i_hat
    model._interpolation_ion_tolerance = epsilon

    model.regenerate_variables()
    model.regenerate_constraints()
    
def apply_ion_weight_constraint(model, i_ref, ion, epsilon):
    # mRNA_ggdw = mRNA_ref
    tot_ion = ion.variable
    mass_coupling_expr = tot_ion - i_ref
    model.add_constraint(kind=InterpolationConstraint,
                        hook=model,
                        id_='ion_interpolation',
                        expr=mass_coupling_expr,
                        lb=-1 * epsilon,
                        ub=epsilon,
                        )
    
def find_rrna_weight_in_rib(model):
    ribosomal_rrna = []
    for rrna in model.rrnas:
        id_ = rrna.id.replace('rrna_','')
        mw_rrna = model.mrnas.get_by_id(id_).molecular_weight
        rrna_w = mw_rrna * symbol_sum([rib.scaling_factor * rib.variable for rib in rrna.ribosomes])
        ribosomal_rrna.append(rrna_w)
    
    return symbol_sum(ribosomal_rrna)
    
def fix_vector_ratio(model, vector, dna_ratio, ppi='ppi_c'):
    '''
    A function similar  to fix_DNA_ratio. It adds a DNA species for the vector
    to the model that with a constant concentration, but this can be used for
    RNAP allocation constraints (to be compatible with those constraints).
    '''
        
    gc_ratio = vector.gc_ratio
    vector_len = vector.len
    vector_dna = DNA(kdeg=0, dna_len=vector_len, gc_ratio=gc_ratio, id="vector")
    # Assumption: kdeg for DNA is close to 0
    vector.add_dna(vector_dna)          
    
    model.add_constraint(kind = ConstantAllocation, 
                                 hook = model, 
                                 expr = vector_dna.variable,
                                 id_ = VECTOR_CONSTANT_CONS_ID,
                                 lb = dna_ratio,
                                 ub = dna_ratio)
    
    # Create dummy DNA reaction
    vector_formation = DNAFormation(id=VECTOR_FORMATION_RXN_ID, name='Vector Formation',
                                 dna=vector_dna, mu_sigma= model._mu_range[-1],
                                 scaled=True)
    model.add_reactions([vector_formation])

    mets = get_dna_synthesis_mets(model, vector_len, gc_ratio, ppi)

    vector_formation.add_metabolites(mets)

    # Add mass balance : 0 = v_syn - [mu]*[DNA]
    model.add_mass_balance_constraint(
        synthesis_flux=vector_formation,
        macromolecule=vector_dna)

    model.regenerate_variables()
    model.regenerate_constraints()