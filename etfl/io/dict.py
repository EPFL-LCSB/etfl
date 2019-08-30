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

import cobra.io.dict as cbd
from cobra.exceptions import SolverNotFound
from optlang.util import expr_to_json, parse_expr
from pytfa.io.dict import get_solver_string, var_to_dict, cons_to_dict, \
    obj_to_dict, rebuild_obj_from_dict
from pytfa.thermo.tmodel import ThermoModel

from ..core.enzyme import Enzyme, Ribosome, Peptide, RNAPolymerase
from ..core.memodel import MEModel
from ..core.rna import mRNA, rRNA, tRNA
from ..core.expression import get_trna_charging_id
from ..core.reactions import TranslationReaction, TranscriptionReaction, \
    EnzymaticReaction, ProteinComplexation, DegradationReaction, ExpressionReaction
from ..core.thermomemodel import ThermoMEModel
from ..optim.utils import rebuild_constraint, rebuild_variable
from ..optim.variables import tRNAVariable
from ..utils.utils import replace_by_enzymatic_reaction, \
    replace_by_translation_reaction, replace_by_transcription_reaction, \
    replace_by_reaction_subclass, replace_by_me_gene

SOLVER_DICT = {
    'optlang.gurobi_interface':'optlang-gurobi',
    'optlang.cplex_interface':'optlang-cplex',
    'optlang.glpk_interface':'optlang-glpk',
}

def metabolite_thermo_to_dict(metthermo):
    return metthermo.thermo.__dict__

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
    return obj


def mrna_to_dict(mrna):
    obj = OrderedDict()
    obj['id'] = mrna.id
    obj['gene_id'] = mrna._gene_id
    obj['kdeg'] = mrna.kdeg
    obj['varname'] = mrna.variable.name
    return obj

def ribosome_to_dict(ribosome):
    obj = enzyme_to_dict(ribosome)
    obj['kribo'] = ribosome.kribo
    obj['rrna_composition'] = ribosome.rrna_composition
    return obj


def rnap_to_dict(rnap):
    obj = enzyme_to_dict(rnap)
    obj['ktrans'] = rnap.ktrans
    obj['kdeg'] = rnap.kdeg
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

def get_solver_string(model):
    return SOLVER_DICT[model.solver.__class__.__module__]


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
            obj['dna_nucleotides'] = model.dna_nucleotides
        except AttributeError:
            # DNA has not been added
            pass
        obj['aa_dict'] = model.aa_dict

        # Growth
        obj['growth_reaction'] = model.growth_reaction.id

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


        obj['kind'] = 'MEModel'
        is_me = True

    if isinstance(model, ThermoMEModel):
        obj['kind'] = 'ThermoMEModel'
        is_me = True


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

        is_peptide = False

        if is_me:
            if the_met_id in model.peptides:
                met_dict['kind'] = 'Peptide'
                met_dict['gene_id'] = the_met._gene_id
                is_peptide = True
            if the_met_id in model.rrnas:
                met_dict['kind'] = 'rRNA'
                is_peptide = True

        if is_thermo and not is_peptide: # peptides have no thermo
            _add_thermo_metabolite_info(the_met, rxn_dict)

    for gene_dict in obj['genes']:
        try:
            gene_dict['sequence'] = str(model.genes.get_by_id(gene_dict['id']).sequence)
        except AttributeError:
            # Not an ExpressedGene
            pass

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

    new._push_queue()

    # Force update GPR info
    for rxn in new.reactions:
        rxn.gene_reaction_rule = rxn.gene_reaction_rule

    for the_var_dict in tqdm(obj['variables'], desc='rebuilding variables'):
        this_id = the_var_dict['id']
        classname = the_var_dict['kind']
        lb = the_var_dict['lb']
        ub = the_var_dict['ub']

        rebuild_variable(classname, new, this_id, lb, ub)

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


    try:
        rebuild_obj_from_dict(new, obj['objective'])
    except KeyError:
        pass

    new.repair()
    return new


def init_me_model_from_dict(new, obj):

    # Convenience attributes
    # new._mu = new.variables.get(obj['_mu'])
    # new.compositions = rebuild_compositions(new, obj['compositions'])
    new.mu_bins = obj['mu_bins']
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
    # new.coupling_dict = rebuild_coupling_dict(new, obj['coupling_dict'])
    new.add_enzymes([enzyme_from_dict(x) for x in obj['enzymes']])

    # Make RNAP
    new_rnap = rnap_from_dict(obj['rnap'])
    new_rnap._model = new
    new.enzymes._replace_on_id(new_rnap)
    new.rnap = new_rnap

    # Make ribosome
    new_rib = ribosome_from_dict(obj['ribosome'])
    new_rib._model = new
    new.enzymes._replace_on_id(new_rib)
    new.ribosome = new_rib

    # Populate Peptides
    find_peptides_from_dict(new, obj)
    find_rrna_from_dict(new, obj)
    # Recover the gene sequences
    find_genes_from_dict(new, obj)


    # Populate mRNAs
    new.add_mrnas([mrna_from_dict(x) for x in obj['mrnas']], add_degradation=False)

    # Populate EnzymaticReaction and TranslationReaction
    find_enzymatic_reactions_from_dict(new, obj)
    find_translation_reactions_from_dict(new, obj)
    find_transcription_reactions_from_dict(new, obj)
    find_complexation_reactions_from_dict(new, obj)
    # link_enzyme_complexation(new, obj)

    # recover tRNAs
    try:
        rebuild_trna(new, obj)
    except KeyError:
        pass

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
    new = init_me_model_from_dict(new,obj)
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
    return mRNA(id = obj['id'],
                  kdeg = obj['kdeg'],
                  gene_id = obj['gene_id'])

def ribosome_from_dict(obj):
    return  Ribosome(   id = obj['id'],
                        kribo = obj['kribo'],
                        kdeg = obj['kdeg'],
                        composition=obj['composition'],
                        rrna=obj['rrna_composition'])

def rnap_from_dict(obj):
    return  RNAPolymerase(  id = obj['id'],
                            ktrans = obj['ktrans'],
                            kdeg = obj['kdeg'],
                            composition=obj['composition'])

def find_enzymatic_reactions_from_dict(new, obj):
    for rxn_dict in obj['reactions']:
        if rxn_dict['kind'] == 'EnzymaticReaction':
            if 'scaled' in rxn_dict:
                scaled = rxn_dict['scaled']
            else:
                scaled = False
            enzymes = [new.enzymes.get_by_id(x) for x in rxn_dict['enzymes']]
            replace_by_enzymatic_reaction(new, rxn_dict['id'],
                                          enzymes,
                                          scaled=scaled)


def find_translation_reactions_from_dict(new, obj):
    new_transl_rxns = list()
    for rxn_dict in obj['reactions']:
        if rxn_dict['kind'] == 'TranslationReaction':
            if 'scaled' in rxn_dict:
                scaled = rxn_dict['scaled']
            else:
                scaled = False
            enzymes = new.ribosome
            enz_rxn = replace_by_translation_reaction(new,
                                                      reaction_id=rxn_dict['id'],
                                                      gene_id=rxn_dict['gene_id'],
                                                      scaled=scaled,
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
            enzymes = new.rnap
            enz_rxn = replace_by_transcription_reaction(new,
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


def find_peptides_from_dict(new, obj):
    new_peptides = list()
    for met_dict in obj['metabolites']:
        if met_dict['kind'] == 'Peptide':
            met = new.metabolites.get_by_id(met_dict['id'])
            pep = Peptide.from_metabolite(met, met_dict['gene_id'])
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
        trna_var_prefix = tRNAVariable.prefix
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

def find_genes_from_dict(new, obj):
    for gene_dict in obj['genes']:
        try:
            sequence = gene_dict['sequence']
            replace_by_me_gene(new, gene_dict['id'], str(sequence))
        except KeyError:
            pass

