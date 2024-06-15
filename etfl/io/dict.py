# -*- coding: utf-8 -*-
"""
.. module:: etfl
   :platform: Unix, Windows
   :synopsis: Expression and thermodynamics-based models

.. moduleauthor:: Pierre Salvy

Make the model serializable
"""
from collections import OrderedDict, defaultdict

from tqdm import tqdm

from Bio.Seq import Seq
#from Bio.Alphabet import DNAAlphabet, RNAAlphabet, ProteinAlphabet

import cobra.io.dict as cbd
from cobra.exceptions import SolverNotFound
from optlang.util import expr_to_json, parse_expr
from pytfa.io.dict import get_solver_string, var_to_dict, cons_to_dict, \
    obj_to_dict, rebuild_obj_from_dict
from pytfa.thermo.tmodel import ThermoModel

from ..core.enzyme import Enzyme, Ribosome, Peptide, RNAPolymerase
from ..core.memodel import MEModel
from ..core.rna import mRNA, rRNA, tRNA
from ..core.dna import DNA
from ..core.genes import ExpressedGene, CodingGene
from ..core.expression import get_trna_charging_id
from ..core.reactions import TranslationReaction, TranscriptionReaction, \
    EnzymaticReaction, ProteinComplexation, DegradationReaction, \
    ExpressionReaction, DNAFormation
from ..core.thermomemodel import ThermoMEModel
from ..optim.utils import rebuild_constraint, rebuild_variable
from ..optim.variables import tRNAVariable, GrowthRate, FreeEnzyme
from ..utils.utils import replace_by_reaction_subclass, replace_by_me_gene, replace_by_coding_gene
from ..core.recombmodel import RecombModel
from ..core.vector import Plasmid

SOLVER_DICT = {
    'optlang.gurobi_interface':'optlang-gurobi',
    'optlang.cplex_interface':'optlang-cplex',
    'optlang.glpk_interface':'optlang-glpk',
}

MW_OVERRIDE_KEY = 'molecular_weight_override'

def metabolite_thermo_to_dict(metthermo):
    return metthermo.thermo.__dict__

def expressed_gene_to_dict(gene):
    obj = OrderedDict()
    obj['id'] = gene.id
    obj['sequence'] = str(gene.sequence)
    obj['rna'] = str(gene.rna)
    obj['copy_number'] = str(gene.copy_number)
    obj['min_tcpt_act'] = str(gene.min_tcpt_activity)
    try:
        obj['transcribed_by'] = [x.id for x in gene.transcribed_by] \
            if gene.transcribed_by is not None else None
    except AttributeError:
        obj['transcribed_by'] = [x for x in gene.transcribed_by] \
            if gene.transcribed_by is not None else None           

    return obj

def coding_gene_to_dict(gene):
    obj = OrderedDict()
    obj['id'] = gene.id
    obj['sequence'] = str(gene.sequence)
    obj['rna'] = str(gene.rna)
    obj['peptide'] = str(gene.peptide)
    obj['copy_number'] = str(gene.copy_number)
    obj['min_tcpt_act'] = str(gene.min_tcpt_activity)
    obj['min_tnsl_act'] = str(gene.min_tnsl_activity)
    try:
        obj['transcribed_by'] = [x.id for x in gene.transcribed_by] \
            if gene.transcribed_by is not None else None
    except AttributeError:
        obj['transcribed_by'] = [x for x in gene.transcribed_by] \
            if gene.transcribed_by is not None else None
    try:
        obj['translated_by'] = [x.id for x in gene.translated_by]\
            if gene.translated_by is not None else None
    except AttributeError:
        obj['translated_by'] = [x for x in gene.translated_by] \
            if gene.translated_by is not None else None            

    return obj

def enzyme_to_dict(enzyme):
    obj = OrderedDict()
    obj['id'] = enzyme.id
    obj['name'] = enzyme.name
    obj['kcat_fwd'] = enzyme.kcat_fwd
    obj['kcat_bwd'] = enzyme.kcat_bwd
    obj['kdeg'] = enzyme.kdeg
    obj['varname'] = enzyme.variable.name
    obj['complexation'] = enzyme.complexation.id
    obj['composition'] = enzyme.composition
    obj['reactions'] = [r.id for r in enzyme.reactions]
    return obj


def mrna_to_dict(mrna):
    obj = OrderedDict()
    obj['id'] = mrna.id
    obj['gene_id'] = mrna._gene_id
    obj['kdeg'] = mrna.kdeg
    obj['varname'] = mrna.variable.name
    if mrna._molecular_weight_override:
        obj['molecular_weight_override'] = mrna._molecular_weight_override
    return obj

def ribosome_to_dict(ribosome):
    if isinstance(ribosome, dict):
        # New models, possibility of several ribosomes
        return {k:_single_ribosome_to_dict(this_ribosome)
                for k,this_ribosome in ribosome.items()}
    else:
        # Older models, there is only one ribosome
        return _single_ribosome_to_dict(ribosome)

def _single_ribosome_to_dict(ribosome):
    obj = enzyme_to_dict(ribosome)
    obj['kribo'] = ribosome.kribo
    obj['rrna_composition'] = ribosome.rrna_composition
    return obj

def rnap_to_dict(rnap):
    if isinstance(rnap, dict):
        # New models, possibility of several rnap
        return {k:_single_rnap_to_dict(this_rnap)
                for k,this_rnap in rnap.items()}
    else:
        # Older models, there is only one rnap
        return _single_rnap_to_dict(rnap)

def _single_rnap_to_dict(rnap):
    obj = enzyme_to_dict(rnap)
    obj['ktrans'] = rnap.ktrans
    obj['kdeg'] = rnap.kdeg
    return obj

def dna_to_dict(dna):
    obj = OrderedDict()
    obj['id'] = dna.id
    obj['name'] = dna.name
    obj['gc_ratio'] = dna.gc_ratio
    obj['len'] = dna.len
    obj['kdeg'] = dna.kdeg
    return obj

def archive_variables(var_dict):
    obj = OrderedDict()

    obj['variables'] = []
    for classname,variables in var_dict.items():
        obj[classname] = list(map(var_to_dict, variables))

    return obj

def archive_constraints(cons_dict):
    obj = OrderedDict()

    for classname,constraints in cons_dict.items():
        obj[classname] = list(map(cons_to_dict, constraints))

    return obj

def archive_compositions(compositions):
    """
    Turns a peptide compositions dict of the form:

    .. code:: python

        { 'b3991': defaultdict(int,
                 {<Metabolite ala__L_c at 0x7f7d25504f28>: -42,
                  <Metabolite arg__L_c at 0x7f7d2550bcf8>: -11,
                  <Metabolite asn__L_c at 0x7f7d2550beb8>: -6,
                  ...}),
        ...}


    to:

    .. code:: python

        { 'b3991': defaultdict(int,,
                {'ala__L_c': -42,
                 'arg__L_c': -11,
                 'asn__L_c': -6,
                  ...}),
        ...}

    :param compositions:
    :return:
    """
    obj = {}
    for k, stoich in compositions.items():
        obj[k] = _stoichiometry_to_dict(stoich)

    return obj

def _stoichiometry_to_dict(stoichiometric_dict):
    """
    Turns a stoichiometric compositions dict of the form:

    .. code:: python

        'b3991': defaultdict(int,
               {<Metabolite ala__L_c at 0x7f7d25504f28>: -42,
                <Metabolite arg__L_c at 0x7f7d2550bcf8>: -11,
                <Metabolite asn__L_c at 0x7f7d2550beb8>: -6,
                ...})

    to:

    .. code:: python

        'b3991': defaultdict(int,,
                {'ala__L_c': -42,
                 'arg__L_c': -11,
                 'asn__L_c': -6,
                  ...})
    """
    return defaultdict(int, {k.id:v for k,v in stoichiometric_dict.items()})

def archive_coupling_dict(coupling_dict):
    """
    Turns an enzyme coupling dict of the form:

    .. code:: python

        {'AB6PGH': <Enzyme AB6PGH at 0x7f7d1371add8>,
         'ABTA': <Enzyme ABTA at 0x7f7d1371ae48>,
         'ACALD': <Enzyme ACALD at 0x7f7d1371aeb8>}

    to:

    .. code:: python

        {'AB6PGH': 'AB6PGH',
         'ABTA': 'ABTA',
         'ACALD': 'ACALD'
    """
    return {k:v.id for k,v in coupling_dict.items()}

def archive_trna_dict(model):
    """
    Turns a tNA information dict of the form:

    .. code:: python

        {'ala__L_c': (<tRNA charged_tRNA_ala__L_c at 0x7f84c16d07b8>,
                      <tRNA uncharged_tRNA_ala__L_c at 0x7f84c16d0be0>,
                      <Reaction trna_ch_ala__L_c at 0x7f84c16d0978>),
         'arg__L_c': (<tRNA charged_tRNA_arg__L_c at 0x7f84c169b588>,
                      <tRNA uncharged_tRNA_arg__L_c at 0x7f84c169b5f8>,
                      <Reaction trna_ch_arg__L_c at 0x7f84c0563ef0>)}

    to:

    .. code:: python

        {'ala__L_c': ('charged_tRNA_ala__L_c',
                      'uncharged_tRNA_ala__L_c',
                      'trna_ch_ala__L_c'),
         'arg__L_c': ('charged_tRNA_arg__L_c',
                      'uncharged_tRNA_arg__L_c',
                      'trna_ch_arg__L_c')}
    """
    return{k:(v[0].id,v[1].id,v[2].id) for k,v in model.trna_dict.items()}

def vector_to_dict(vector_dict, vector_copy_number):
    if isinstance(vector_dict, dict):
        # New models, possibility of several ribosomes
        return {k:_single_vector_to_dict(this_vector, vector_copy_number[k])
                for k,this_vector in vector_dict.items()}
    else:
        # Older models, there is only one ribosome
        return _single_vector_to_dict(vector_dict)
    
def _single_vector_to_dict(vector, copy_number):
    obj = OrderedDict()
    obj['id'] = vector.id
    obj['genes'] = [gene.id for gene in vector.genes]
    obj['reactions'] = [rxn.id for rxn in vector.reactions]
    obj['proteins'] = [prot.id for prot in vector.proteins]
    obj['peptides'] = [pep.id for pep in vector.peptides]
    obj['rnap'] = [x.id for x in vector.rnap] \
            if vector.rnap is not None else None
    obj['ribosome'] = [x.id for x in vector.ribosome] \
            if vector.ribosome is not None else None
    obj['gc_ratio'] = vector.gc_ratio
    obj['len'] = vector.len
    obj['sequence'] = vector.sequence if vector.sequence is not None else None
    obj['default_rnap'] = vector.default_rnap
    obj['default_rib'] = vector.default_rib
    obj['calibration_tcpt'] = vector.calibration_tcpt
    obj['calibration_tnsl'] = vector.calibration_tnsl
    obj['copy_number'] = copy_number
    return obj
    

def model_to_dict(model):
    """

    :param model:
    :return:
    """

    # Take advantage of cobra's dict serialization for metabolites and
    # reactions
    obj = cbd.model_to_dict(model)

    obj['solver'] = get_solver_string(model)
    obj['objective'] = obj_to_dict(model)

    # Copy variables, constraints
    # obj['var_dict'] = archive_variables(model._var_kinds)
    # obj['cons_dict'] = archive_constraints(model._cons_kinds)
    obj['variables'] = list(map(var_to_dict, model._var_dict.values()))
    obj['constraints'] = list(map(cons_to_dict, model._cons_dict.values()))

    is_thermo = False
    is_me = False

    if isinstance(model, ThermoModel):
        obj['kind'] = 'ThermoModel'
        obj['thermo_data'] = model.thermo_data #it's a dict
        obj['name'] = model.name
        obj['temperature'] = model.TEMPERATURE
        obj['min_ph'] = model.MIN_pH
        obj['max_ph'] = model.MAX_pH
        is_thermo = True

        # Relaxation info
        try:
            obj['relaxation'] = model.relaxation
        except AttributeError:
            pass

    if isinstance(model, MEModel):

        # Convenience attributes
        # obj['_mu'] = model.mu.name
        # obj['compositions'] = archive_compositions(model.compositions)
        # obj['coupling_dict'] = archive_coupling_dict(model.coupling_dict)
        obj['mu_bins'] = model.mu_bins
        obj['essentials'] = model.essentials
        obj['rna_nucleotides'] = model.rna_nucleotides
        obj['rna_nucleotides_mp'] = model.rna_nucleotides_mp
        try:
            obj['dna_nucleotides'] = model.dna_nucleotides
        except AttributeError:
            # DNA has not been added
            pass
        obj['aa_dict'] = model.aa_dict

        # Growth
        obj['growth_reaction'] = model.growth_reaction.id

        # Genes
        obj['expressed_genes'] = list(map(expressed_gene_to_dict,
                                [g for g in model.genes
                                 if type(g) == ExpressedGene]))
        obj['coding_genes'] = list(map(coding_gene_to_dict,
                                [g for g in model.genes
                                 if type(g) == CodingGene]))

        # Enzymes
        obj['enzymes'] = list(map(enzyme_to_dict, model.enzymes))

        # mRNAs
        obj['mrnas'] = list(map(mrna_to_dict, model.mrnas))

        # tRNAs
        obj['trna_dict'] = archive_trna_dict(model)

        # Ribosome
        obj['ribosome'] = ribosome_to_dict(model.ribosome)

        # RNAP
        obj['rnap'] = rnap_to_dict(model.rnap)

        try:
            obj['dna'] = dna_to_dict(model.dna)
        except AttributeError:
            # DNA has not been added
            pass

        obj['kind'] = 'MEModel'
        is_me = True

    if isinstance(model, ThermoMEModel):
        obj['kind'] = 'ThermoMEModel'
        is_me = True
        
    if isinstance(model, RecombModel):
        obj['vector'] = vector_to_dict(model.vectors, model.vector_copy_number)
        obj['biomass_metabolites'] = model.biomass_metabolites
        obj['biomass_composition'] = model.biomass_composition
        obj['kind'] = 'RecombModel'
        

    # Metabolite and Reaction-level cleanup
    for rxn_dict in obj['reactions']:
        rxn = model.reactions.get_by_id(rxn_dict['id'])

        if is_me:
            _add_me_reaction_info(rxn, rxn_dict)

        if is_thermo:
            _add_thermo_reaction_info(rxn, rxn_dict)


    # Peptides and Thermo
    for met_dict in obj['metabolites']:
        the_met_id = met_dict['id']

        met_dict['kind'] = 'Metabolite'
        the_met = model.metabolites.get_by_id(the_met_id)

        is_me_met = False

        if is_me:
            if the_met_id in model.peptides:
                met_dict['kind'] = 'Peptide'
                met_dict['gene_id'] = the_met._gene_id
                if the_met._molecular_weight_override:
                    met_dict['molecular_weight_override'] = \
                        the_met._molecular_weight_override
                is_me_met = True
            if the_met_id in model.rrnas:
                met_dict['kind'] = 'rRNA'
                is_me_met = True

        if is_thermo and not is_me_met: # peptides have no thermo
            _add_thermo_metabolite_info(the_met, rxn_dict)

    return obj

def _add_me_reaction_info(rxn, rxn_dict):

    if isinstance(rxn, ExpressionReaction):
        rxn_dict['scaled'] = rxn._scaled
    # We start with translation reactions because they are also
    # enzymatic reactions
    # Translation Reactions
    if isinstance(rxn, TranslationReaction):
        rxn_dict['kind'] = 'TranslationReaction'
        rxn_dict['gene_id'] = rxn.gene.id
        rxn_dict['trna_stoich'] = rxn.trna_stoich
    # Transcription Reactions
    elif isinstance(rxn, TranscriptionReaction):
        rxn_dict['kind'] = 'TranscriptionReaction'
        rxn_dict['gene_id'] = rxn.gene.id
    # Protein Complexation
    elif isinstance(rxn, ProteinComplexation):
        rxn_dict['kind'] = 'ProteinComplexation'
        rxn_dict['gene_id'] = None
        rxn_dict['target'] = rxn.target.id
    # Degradation Reaction
    elif isinstance(rxn, DegradationReaction):
        rxn_dict['kind'] = 'DegradationReaction'
        rxn_dict['gene_id'] = None
        rxn_dict['macromolecule'] = rxn.macromolecule.id
        rxn_dict['macromolecule_kind'] = rxn.macromolecule.__class__.__name__
    # Enzymatic Reactions
    elif isinstance(rxn, EnzymaticReaction):
        rxn_dict['kind'] = 'EnzymaticReaction'
        rxn_dict['enzymes'] = [x.id for x in rxn.enzymes]
    elif isinstance(rxn, DNAFormation):
        rxn_dict['kind'] = 'DNAFormation'
    # Generic Reaction
    else:
        rxn_dict['kind'] = 'Reaction'

def _add_thermo_reaction_info(rxn, rxn_dict):
    if hasattr(rxn, 'thermo'):
        rxn_dict['thermo'] = rxn.thermo

def _add_thermo_metabolite_info(met, met_dict):
    if hasattr(met, 'thermo'):
        met_dict['thermo'] = metabolite_thermo_to_dict(met)

def model_from_dict(obj, solver=None):
    # Take advantage of cobra's serialization of mets and reactions
    new = cbd.model_from_dict(obj)

    if solver is not None:
        try:
            new.solver = solver
        except SolverNotFound as snf:
            raise snf
    else:
        try:
            new.solver = obj['solver']
        except KeyError:
            pass

    # Populate variables and constraints
    if obj['kind'] == 'ThermoMEModel':
        new = ThermoMEModel(thermo_data=obj['thermo_data'],
                            model=new,
                            name=obj['name'],
                            temperature=obj['temperature'],
                            min_ph=obj['min_ph'],
                            max_ph=obj['max_ph'])
        new = init_thermo_me_model_from_dict(new, obj)
    elif obj['kind'] == 'MEModel':
        # Cast to MEModel
        new = MEModel(new)
        new = init_me_model_from_dict(new, obj)
    elif obj['kind'] == 'ThermoModel':
        new = ThermoModel(thermo_data=obj['thermo_data'],
                          model=new,
                          name=obj['name'],
                          temperature=obj['temperature'],
                          min_ph=obj['min_ph'],
                          max_ph=obj['max_ph'])
        new = init_thermo_model_from_dict(new, obj)
    elif obj['kind'] == 'RecombModel':
        new = MEModel(new)
        new = init_me_model_from_dict(new, obj)
        new = RecombModel(new, inplace=True)
        new = init_recomb_model_from_dict(new, obj)
        

    new._push_queue()

    # Force update GPR info
    for rxn in new.reactions:
        rxn.gene_reaction_rule = rxn.gene_reaction_rule

    for the_var_dict in tqdm(obj['variables'], desc='rebuilding variables'):
        this_id = the_var_dict['id']
        classname = the_var_dict['kind']
        lb = the_var_dict['lb']
        ub = the_var_dict['ub']
        try: #Backward compat
            scaling_factor = the_var_dict['scaling_factor']
        except KeyError:
            scaling_factor = 1

        rebuild_variable(classname, new, this_id, lb, ub, scaling_factor)

    new._push_queue()

    variable_parse_dict = {x.name:x for x in new.variables}

    for the_cons_dict in tqdm(obj['constraints'], desc='rebuilding constraints'):
        this_id = the_cons_dict['id']
        classname = the_cons_dict['kind']
        new_expr = parse_expr(the_cons_dict['expression'],
                              local_dict = variable_parse_dict)
        # str_expr = the_cons_dict['expression']
        #
        # # Sympify the expression so that we can substitute variables afterwards
        # sym_expr = sympify(str_expr)
        #
        # subs_dict = {x:new.variables.get(x.name) for x in sym_expr.free_symbols}
        #
        # new_expr = sym_expr.subs(subs_dict)
        lb = the_cons_dict['lb']
        ub = the_cons_dict['ub']

        rebuild_constraint(classname, new, this_id, new_expr, lb, ub)


    # Mu variable handle for ME-models
    if obj['kind'] in ['ThermoMEModel','MEModel','RecombModel']:
        prostprocess_me(new)

    try:
        rebuild_obj_from_dict(new, obj['objective'])
    except KeyError:
        pass

    new.update_enzyme_reaction_rel()
    new.repair()
    return new


def prostprocess_me(new):
    new._mu = new.get_variables_of_type(GrowthRate).get_by_id('total')


def init_me_model_from_dict(new, obj):

    # Convenience attributes
    # new._mu = new.variables.get(obj['_mu'])
    # new.compositions = rebuild_compositions(new, obj['compositions'])
    new.mu_bins = obj['mu_bins']
    new._mu_range = [new.mu_bins[0][1][0], new.mu_bins[-1][1][-1]]
    new.essentials = obj['essentials']
    new.rna_nucleotides = obj['rna_nucleotides']
    new.rna_nucleotides_mp = obj['rna_nucleotides_mp']
    try:
        new.dna_nucleotides = obj['dna_nucleotides']
    except KeyError:
        # No DNA data
        pass
    new.aa_dict = obj['aa_dict']
    # new.trna_dict = obj['trna_dict']

    # Add growth reaction
    new.growth_reaction = obj['growth_reaction']

    # Populate enzymes
    # new.coupling_dict = rebuild_coupling_dict(, obj['coupling_dict'])
    new.add_enzymes([enzyme_from_dict(x) for x in obj['enzymes']], prep=False)

    # Recover the gene sequences
    find_genes_from_dict(new, obj)
    
    # Make RNAP
    new_rnap = rnap_from_dict(obj['rnap'])
    for this_rnap in new_rnap.values():
        this_rnap._model = new
        new.enzymes._replace_on_id(this_rnap)
    new.rnap = new_rnap

    # Make ribosome
    new_rib = ribosome_from_dict(obj['ribosome'])
    for this_rib in new_rib.values():
        this_rib._model = new
        new.enzymes._replace_on_id(this_rib)
        # adding rrnas to the model
        new.add_rrnas_to_rib_assembly(this_rib)
    new.ribosome = new_rib

    try:
        new.add_dna(dna_from_dict(obj['dna']))
    except KeyError:
        # There is no DNA in the model
        pass

    # Populate Peptides
    find_peptides_from_dict(new, obj)
    find_rrna_from_dict(new, obj)
    

    # Populate mRNAs
    new.add_mrnas([mrna_from_dict(x) for x in obj['mrnas']], add_degradation=False)

    # recover tRNAs
    try:
        rebuild_trna(new, obj)
    except KeyError:
        pass

    # Populate EnzymaticReaction and TranslationReaction
    find_enzymatic_reactions_from_dict(new, obj)
    find_translation_reactions_from_dict(new, obj)
    find_transcription_reactions_from_dict(new, obj)
    find_complexation_reactions_from_dict(new, obj)
    find_dna_formation_reaction_from_dict(new, obj)
    # link_enzyme_complexation(new, obj)

    # Finally, add degradations
    find_degradation_reactions_from_dict(new, obj)

    # new.rnap.init_variable()
    # new.ribosome.init_variable()
    for enz in new.enzymes:
        enz.init_variable()
    for mrna in new.mrnas:
        mrna.init_variable()
    # This is already done by model.add_trna
    # for trna in new.trnas:
    #     trna.init_variable()
    # new.init_ribosome_variables()

    return new

def init_thermo_model_from_dict(new, obj):
    for rxn_dict in obj['reactions']:
        rxn = new.reactions.get_by_id(rxn_dict['id'])

        if 'thermo' in rxn_dict:
            rxn.thermo = rxn_dict['thermo']

    for met_dict in obj['metabolites']:
        met = new.metabolites.get_by_id(met_dict['id'])

        if 'thermo' in met_dict:
            new._prepare_metabolite(met)

    # Relaxation info
    try:
        new.relaxation = obj['relaxation']
    except KeyError:
        pass

    return new


def init_thermo_me_model_from_dict(new, obj):
    new = init_thermo_model_from_dict(new, obj)
    new = init_me_model_from_dict(new, obj)
    return new

def init_recomb_model_from_dict(new, obj):
    new.biomass_metabolites = obj['biomass_metabolites']
    new.biomass_composition = obj['biomass_composition']
    find_vector_from_dict(new, obj['vector'])
    return new

def rebuild_compositions(new, compositions_dict):
    """
    Performs the reverse operation of :func:archive_compositions

    :param new:
    :param compositions_dict:
    :return:
    """

    return {k:_rebuild_stoichiometry(new,stoich)
            for k,stoich in compositions_dict.items()}

def _rebuild_stoichiometry(new, stoich):
    """
    Performs the reverse operation of :func:_stoichiometry_to_dict

    :param new:
    :param stoich:
    :return:
    """

    return defaultdict(int,
                       {new.metabolites.get_by_id(k):v
                        for k,v in stoich.items()})

def rebuild_coupling_dict(new, coupling_dict):
    """
    Performs the reverse operation of :func:archive_coupling_dict

    :param new:
    :param coupling_dict:
    :return:
    """
    return {k:new.enzymes.get_by_id(v) for k,v in coupling_dict.items()}

def enzyme_from_dict(obj):
    return Enzyme(id = obj['id'],
                     kcat_fwd = obj['kcat_fwd'],
                     kcat_bwd = obj['kcat_bwd'],
                     kdeg = obj['kdeg'],
                     composition=obj['composition'])

def mrna_from_dict(obj):
    the_mrna = mRNA(id = obj['id'],
                    kdeg = obj['kdeg'],
                    gene_id = obj['gene_id'])
    if MW_OVERRIDE_KEY in obj:
        the_mrna._molecular_weight_override = obj[MW_OVERRIDE_KEY]
    return the_mrna

def ribosome_from_dict(obj):
    if isinstance(list(obj.values())[0],dict):
        # The object describes a collection of ribosomes, it is one of the newer
        # models which allow several ribosomes
        return {k : _single_ribosome_from_dict(x) for k,x in obj.items()}
    else:
        # This is an older model with a single ribosome
        return {obj['id']:_single_ribosome_from_dict(obj)}

def _single_ribosome_from_dict(obj):
    return  Ribosome(   id = obj['id'],
                        kribo = obj['kribo'],
                        kdeg = obj['kdeg'],
                        composition=obj['composition'],
                        rrna=obj['rrna_composition'])

def rnap_from_dict(obj):
    if isinstance(list(obj.values())[0],dict):
        # The object describes a collection of rnap, it is one of the newer
        # models which allow several rnap
        return {k : _single_rnap_from_dict(x) for k,x in obj.items()}
    else:
        # This is an older model with a single rnap
        return {obj['id']:_single_rnap_from_dict(obj)}

def _single_rnap_from_dict(obj):
    return  RNAPolymerase(  id = obj['id'],
                            ktrans = obj['ktrans'],
                            kdeg = obj['kdeg'],
                            composition=obj['composition'])

def dna_from_dict(obj):
    return DNA( id = obj['id'],
                dna_len = obj['len'],
                name=obj['name'],
                gc_ratio = obj['gc_ratio'],
                kdeg=obj['kdeg'])

def find_enzymatic_reactions_from_dict(new, obj):
    for rxn_dict in obj['reactions']:
        if rxn_dict['kind'] == 'EnzymaticReaction':
            if 'scaled' in rxn_dict:
                scaled = rxn_dict['scaled']
            else:
                scaled = False
            enzymes = [new.enzymes.get_by_id(x) for x in rxn_dict['enzymes']]
            replace_by_reaction_subclass(new,
                                          kind=EnzymaticReaction,
                                          reaction_id=rxn_dict['id'],
                                          enzymes=enzymes,
                                          scaled=scaled)


def find_translation_reactions_from_dict(new, obj):
    new_transl_rxns = list()
    for rxn_dict in obj['reactions']:
        if rxn_dict['kind'] == 'TranslationReaction':
            if 'scaled' in rxn_dict:
                scaled = rxn_dict['scaled']
            else:
                scaled = False
            if 'trna_stoich' in rxn_dict:
                trna_stoich = rxn_dict['trna_stoich']
            else:
                trna_stoich = None

            enzymes = new.ribosome.values()
            enz_rxn = replace_by_reaction_subclass(new,
                                                      kind = TranslationReaction,
                                                      reaction_id=rxn_dict['id'],
                                                      gene_id=rxn_dict['gene_id'],
                                                      scaled=scaled,
                                                      trna_stoich=trna_stoich,
                                                      enzymes=enzymes)
            new_transl_rxns.append(enz_rxn)
    new.translation_reactions += new_transl_rxns


def find_transcription_reactions_from_dict(new, obj):
    new_transc_rxns = list()
    for rxn_dict in obj['reactions']:
        if rxn_dict['kind'] == 'TranscriptionReaction':
            if 'scaled' in rxn_dict:
                scaled = rxn_dict['scaled']
            else:
                scaled = False
            enzymes = new.rnap.values()
            enz_rxn = replace_by_reaction_subclass(new,
                                                   kind=TranscriptionReaction,
                                                        reaction_id=rxn_dict['id'],
                                                        gene_id=rxn_dict['gene_id'],
                                                        scaled=scaled,
                                                        enzymes=enzymes)
            new_transc_rxns.append(enz_rxn)
    new.transcription_reactions += new_transc_rxns


def find_complexation_reactions_from_dict(new, obj):
    new_rxns = list()
    for rxn_dict in obj['reactions']:
        if rxn_dict['kind'] == 'ProteinComplexation':
            target = new.enzymes.get_by_id(rxn_dict['target'])
            if 'scaled' in rxn_dict:
                scaled = rxn_dict['scaled']
            else:
                scaled = False
            new_rxn = replace_by_reaction_subclass(new,
                                                   kind = ProteinComplexation,
                                                   reaction_id=rxn_dict['id'],
                                                   scaled=scaled,
                                                   target=target)
            new_rxns.append(new_rxn)
    new.complexation_reactions += new_rxns

def link_enzyme_complexation(new, obj):
    for enz_obj in obj['enzymes']:
        the_enz = new.enzymes.get_by_id(enz_obj['id'])
        the_enz.complexation = new.reactions.get_by_id(enz_obj['complexation'])

def find_degradation_reactions_from_dict(new, obj):
    new_rxns = list()
    for rxn_dict in obj['reactions']:
        if rxn_dict['kind'] == 'DegradationReaction':
            if rxn_dict['macromolecule_kind']    == 'mRNA':
                macromol = new.mrnas.get_by_id(rxn_dict['macromolecule'])
            elif  rxn_dict['macromolecule_kind'] == 'Enzyme' or \
                  rxn_dict['macromolecule_kind'] == 'RNAPolymerase' or \
                  rxn_dict['macromolecule_kind'] == 'Ribosome':
                macromol = new.enzymes.get_by_id(rxn_dict['macromolecule'])
            else:
                raise(TypeError('Macromolecule type not recognized: {}'
                      .format(rxn_dict['macromolecule_kind'])))
            if 'scaled' in rxn_dict:
                scaled = rxn_dict['scaled']
            else:
                scaled = False
            new_rxn = replace_by_reaction_subclass(new,
                                                   kind = DegradationReaction,
                                                   reaction_id=rxn_dict['id'],
                                                   scaled=scaled,
                                                   macromolecule=macromol)
            new_rxns.append(new_rxn)
    new.degradation_reactions += new_rxns

def find_dna_formation_reaction_from_dict(new, obj):
    new_rxns = list()
    for rxn_dict in obj['reactions']:
        if rxn_dict['kind'] == 'DNAFormation':
            dna = new.dna
            new_rxn = replace_by_reaction_subclass(new,
                                                   kind = DNAFormation,
                                                   reaction_id=rxn_dict['id'],
                                                   scaled=rxn_dict['scaled'],
                                                   dna=dna,
                                                   mu_sigma=new._mu_range[-1])
            new_rxns.append(new_rxn)

def find_peptides_from_dict(new, obj):
    new_peptides = list()
    for met_dict in obj['metabolites']:
        if met_dict['kind'] == 'Peptide':
            met = new.metabolites.get_by_id(met_dict['id'])
            pep = Peptide.from_metabolite(met, met_dict['gene_id'])
            if MW_OVERRIDE_KEY in met_dict:
                pep._molecular_weight_override = met_dict[MW_OVERRIDE_KEY]
            new.metabolites._replace_on_id(pep)
            new_peptides.append(pep)
    new.peptides += new_peptides

def find_rrna_from_dict(new, obj):
    new_rrnas = list()
    for met_dict in obj['metabolites']:
        if met_dict['kind'] == 'rRNA':
            met = new.metabolites.get_by_id(met_dict['id'])
            rrna = rRNA.from_metabolite(met)
            new.metabolites._replace_on_id(rrna)
            new_rrnas.append(rrna)
    new.rrnas += new_rrnas

def rebuild_trna(new, obj):
        new_trna_dict = obj['trna_dict']
        new_trnas = list()
        for aa_id, (ch_trna_id, unch_trna_id, charging_rxn_id) in new_trna_dict.items():
            aa = new.metabolites.get_by_id(aa_id)

            charged_trna = tRNA(aminoacid_id=aa.id,
                                charged=True,
                                name=aa.name)
            charged_trna.id = ch_trna_id

            uncharged_trna = tRNA(aminoacid_id=aa.id,
                                  charged=False,
                                  name=aa.name)
            uncharged_trna.id = unch_trna_id

            new_trnas.extend([charged_trna,uncharged_trna])

            charging_reaction = new.reactions.get_by_id(get_trna_charging_id(aa_id))

            new_trna_dict[aa_id] = (charged_trna, uncharged_trna, charging_reaction)

        new.add_trnas(new_trnas)
        new.trna_dict = new_trna_dict

def find_genes_from_dict(new, obj):

    for gene_dict in obj['coding_genes']:
        try:
            sequence = gene_dict['sequence']
            g = replace_by_me_gene(new, gene_dict['id'], str(sequence))
            g._rna            = Seq(gene_dict['rna'])#, alphabet=RNAAlphabet())
            g._peptide        = Seq(gene_dict['peptide'])#, alphabet=ProteinAlphabet())
            g._copy_number    = int(gene_dict['copy_number'])
            g._transcribed_by = [new.enzymes.get_by_id(e)
                                 for e in gene_dict['transcribed_by']] \
                                 if gene_dict['transcribed_by'] else None
            g._translated_by  = [new.enzymes.get_by_id(e)
                                 for e in gene_dict['translated_by']] \
                                 if gene_dict['translated_by'] else None
        except KeyError:
            pass
        
    for gene_dict in obj['expressed_genes']:
        try:
            sequence = gene_dict['sequence']
            g = replace_by_me_gene(new, gene_dict['id'], str(sequence))
            g = replace_by_coding_gene(new, gene_dict['id'])
            g._rna            = Seq(gene_dict['rna'])#, alphabet=RNAAlphabet())
            g._copy_number    = int(gene_dict['copy_number'])
            g._transcribed_by = [new.enzymes.get_by_id(e)
                                 for e in gene_dict['transcribed_by']] \
                                 if gene_dict['transcribed_by'] else None
        except KeyError:
            pass
        
def find_vector_from_dict(new, obj):
    new.vector_copy_number = dict()
    new.vectors = dict()
    if isinstance(list(obj.values())[0],dict):
        # The object describes a collection of ribosomes, it is one of the newer
        # models which allow several ribosomes
        return {k : _single_vector_from_dict(new, x) for k,x in obj.items()}
    else:
        # This is an older model with a single ribosome
        return {obj['id']:_single_vector_from_dict(new, obj)}
    
def _single_vector_from_dict(new, obj):
    
    vector = Plasmid(id_ = obj['id'], 
                     genes = [new.genes.get_by_id(x) for x in obj['genes']], 
                     reactions = [new.reactions.get_by_id(x) for x in obj['reactions']], 
                     gc_ratio = obj['gc_ratio'], 
                     length = obj['len'],
                     proteins = [new.enzymes.get_by_id(x) for x in obj['proteins']],
                     sequence = obj['sequence'] if obj['sequence'] else None,
                     mrna_dict = None,
                     rnap = [new.enzymes.get_by_id(e) for e in obj['rnap']] \
                         if obj['rnap'] else None,
                     ribosome = [new.enzymes.get_by_id(e) for e in obj['ribosome']] \
                         if obj['ribosome'] else None)
        
    vector.peptides = [new.peptides.get_by_id(p) for p in obj['peptides']]
    vector.default_rnap = new.enzymes.get_by_id(obj['default_rnap'])
    vector.default_rib = new.enzymes.get_by_id(obj['default_rib'])
    vector.calibration_tcpt = obj['calibration_tcpt']
    vector.calibration_tnsl = obj['calibration_tnsl']
    # Just to be added manually to the model
    new.vectors[vector.id] = vector
    # assigning the vector copy number
    new.vector_copy_number[vector.id] = obj['copy_number']
