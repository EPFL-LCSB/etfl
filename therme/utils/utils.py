from therme.therme.core.reaction import EnzymaticReaction, TranslationReaction

def replace_by_enzymatic_reaction(model, reaction_id, enzymes):
    rxn = model.reactions.get_by_id(reaction_id)
    enz_rxn = EnzymaticReaction.from_reaction(reaction=rxn,
                                              enzymes=enzymes)
    _replace_by_me_reaction(model, rxn, enz_rxn)
    return enz_rxn


def replace_by_translation_reaction(model, reaction_id, enzymes):
    rxn = model.reactions.get_by_id(reaction_id)
    enz_rxn = TranslationReaction.from_reaction(reaction=rxn,
                                              enzymes=enzymes)
    _replace_by_me_reaction(model, rxn, enz_rxn)
    return enz_rxn


def _replace_by_me_reaction(model, rxn, enz_rxn):

    # model.remove_reactions(reactions=[rxn])
    # model.add_reactions([enz_rxn])
    model.reactions._replace_on_id(enz_rxn)
    enz_rxn._model = model
    enz_rxn.notes = rxn.notes
    if hasattr(rxn, 'thermo'):
        enz_rxn.thermo = rxn.thermo