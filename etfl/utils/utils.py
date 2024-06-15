from ..core.expression import _extract_trna_from_reaction, make_stoich_from_aa_sequence
from ..core.reactions import EnzymaticReaction, TranslationReaction, TranscriptionReaction
from ..core.genes import ExpressedGene, CodingGene

def replace_by_enzymatic_reaction(model, reaction_id, enzymes, scaled):
    rxn = model.reactions.get_by_id(reaction_id)
    enz_rxn = EnzymaticReaction.from_reaction(reaction=rxn,
                                              enzymes=enzymes,
                                              scaled=scaled)
    _replace_by_me_reaction(model, rxn, enz_rxn)
    return enz_rxn


def replace_by_translation_reaction(model, reaction_id, gene_id, enzymes,
                                    trna_stoich, scaled):
    rxn  = model.reactions.get_by_id(reaction_id)
    enz_rxn = TranslationReaction.from_reaction(reaction=rxn,
                                                gene_id=gene_id,
                                                enzymes=enzymes,
                                                trna_stoich=trna_stoich,
                                                scaled=scaled)
    _replace_by_me_reaction(model, rxn, enz_rxn)
    return enz_rxn

def replace_by_transcription_reaction(model, reaction_id, gene_id, enzymes, scaled):
    rxn  = model.reactions.get_by_id(reaction_id)
    enz_rxn = TranscriptionReaction.from_reaction(reaction=rxn,
                                                  gene_id=gene_id,
                                                  enzymes=enzymes,
                                                  scaled=scaled)
    _replace_by_me_reaction(model, rxn, enz_rxn)
    return enz_rxn

def replace_by_reaction_subclass(model, kind, reaction_id, **kwargs):
    rxn  = model.reactions.get_by_id(reaction_id)
    new_rxn = kind.from_reaction(reaction=rxn, **kwargs)
    _replace_by_me_reaction(model, rxn, new_rxn)
    return new_rxn


def _replace_by_me_reaction(model, rxn, enz_rxn):

    # model.remove_reactions(reactions=[rxn])
    # model.add_reactions([enz_rxn])
    model.reactions._replace_on_id(enz_rxn)
    enz_rxn._model = model

    if hasattr(rxn, 'thermo'):
        enz_rxn.thermo = rxn.thermo


def replace_by_me_gene(model, gene_id, sequence):
    gene = model.genes.get_by_id(gene_id)
    new = CodingGene.from_gene(gene=gene,
                                  sequence=sequence)

    # # That is not a typo, see class cobra.core.Species
    # if gene.reactions:
    #     new._reaction = set(x for x in gene.reactions)
    # else:
    #     # Find by enumerating:
    #     new._reaction = set(r for r in model.reactions
    #                         if gene.id in [x.id for x in r.genes])

    model.genes._replace_on_id(new)
    new._model = model
    new.notes = gene.notes
    return new

def replace_by_coding_gene(model, gene_id):
    # a function to convert coding gene to expressed gene
    gene = model.genes.get_by_id(gene_id)
    new = ExpressedGene(id=gene_id , name=gene_id , 
                     sequence=gene.sequence)


    model.genes._replace_on_id(new)
    new._model = model
    new.notes = gene.notes
    return new
