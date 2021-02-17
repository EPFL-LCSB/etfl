import sympy

from ..core.enzyme import Enzyme
from ..utils.parsing import parse_gpr

def infer_enzyme_from_gpr(reaction, default_kcat, default_kdeg):
    new_enzymes = list()
    compositions = compositions_from_gpr(reaction)
    for e,composition in enumerate(compositions):
        new_enzyme = Enzyme(id = reaction.id + '_inferred_{}'.format(e),
                            kcat=default_kcat,
                            kdeg=default_kdeg,
                            composition=composition)
        new_enzymes.append(new_enzyme)
    return new_enzymes

def compositions_from_gpr(reaction):
    """
    *Warning*: Use this function only if you have no information on the enzymes.
    Logically parses the GPR to automatically find isozymes ( logical OR )
    and subunits ( logical AND ), and creates the necessary complexation
    reactions: 1 per isozyme, requiring the peptides of each subunit

    :param reaction:
    :type reaction: cobra.Reaction
    :return:
    """

    model = reaction.model

    this_gpr = reaction.gene_reaction_rule

    sym_gpr = parse_gpr(this_gpr)

    if isinstance(sym_gpr, sympy.Symbol):
        # GPR of the type: '(gene0)'
        # Gene <=> Protein
        isozymes = [sym_gpr]
    elif isinstance(sym_gpr, sympy.And):
        # GPR of the type: '(gene0 & gene1)'
        # Subunits of one enzyme
        isozymes = [sym_gpr]
    elif isinstance(sym_gpr, sympy.Or):
        # GPR of the type: '(gene0 | gene1)', '((gene0 & gene1) | gene2)'
        # Two isozymes that are the arguments of the OR
        isozymes = sym_gpr.args

    compositions = []

    for e, this_isozyme in enumerate(isozymes):
        if isinstance(this_isozyme, sympy.And):
            # this is a GPR with several subunits
            peptides = {x.name: 1 \
                        for x in this_isozyme.args}
        elif isinstance(this_isozyme, sympy.Symbol):
            # there is only one subunit
            peptides = {this_isozyme.name: 1}
        else:
            # The GPR has been incorrectly parsed
            model.logger.error('Incorrect parsing of {}'.format(isozymes))
            raise TypeError

        compositions += [peptides]

    return compositions