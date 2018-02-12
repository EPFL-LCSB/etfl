# -*- coding: utf-8 -*-
"""
.. module:: thermome
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

ME-related Reaction subclasses and methods definition


"""
from cobra import Reaction, DictList

class EnzymaticReaction(Reaction):
    """
    Subclass to describe reactions that are catalyzed by an enzyme.
    """
    def __init__(self, enzymes = None, *args, **kwargs):
        Reaction.__init__(self, *args, **kwargs)
        self.enzymes = DictList()
        if enzymes:
            self.add_enzymes(enzymes)

    @staticmethod
    def from_reaction(reaction, enzymes = None):
        """
        This method clones a cobra.Reaction object into an enzymatic reaction,
        and attaches enzymes to it

        :param reaction: the reaction to reproduce
        :type reaction: cobra.Reaction
        :param enzymes: a(n iterable of the) enzyme(s) to be attached to the reaction
        :return: an EnzymaticReaction object
        """
        new =  EnzymaticReaction(id = reaction.id,
                                 name= reaction.name,
                                 subsystem= reaction.subsystem,
                                 lower_bound= reaction.lower_bound,
                                 upper_bound= reaction.upper_bound,
                                 enzymes= enzymes)
        new.add_metabolites(reaction.metabolites)
        new.gene_reaction_rule = reaction.gene_reaction_rule
        return new


    def add_enzymes(self, enzymes):
        """
        Method to add the enzymes to the reaction.
        :param enzymes: iterable of or single Enzyme object
        :return:
        """

        if not hasattr(enzymes, '__iter__'):
            enzymes = [enzymes]

        for e in enzymes:
            # Avoid duplicates
            if e.id in self.enzymes:
                self.enzymes._replace_on_id(e)
            else:
                self.enzymes += [e]


class TranscriptionReaction(EnzymaticReaction):
    """
    Class describing translation - Assembly of amino acids into peptides
    """

    def __init__(self, id, name, sequence, enzyme, **kwargs):
        EnzymaticReaction.__init__(self,
                                   id=id,
                                   name=name,
                                   enzymes=enzyme,
                                   **kwargs)
        self.sequence = sequence

    # The number of amino acids is in the stoichiometry
    @property
    def nucleotide_length(self):
        return len(self.sequence)


    def add_rnap(self, ribosome):
        """
        By definition this reaction will be catalyzed by RNA polymerase
        :param ribosome:
        :type ribosome: pytfa.me.Ribosome
        :return:
        """
        self.enzyme = ribosome

class TranslationReaction(EnzymaticReaction):
    """
    Class describing translation - Assembly of amino acids into peptides
    """

    def __init__(self, id, name, sequence, enzyme, **kwargs):
        EnzymaticReaction.__init__(self,
                                   id=id,
                                   name=name,
                                   enzymes=enzyme,
                                   **kwargs)
        self.sequence = sequence

    # The number of amino acids is in the stoichiometry
    @property
    def aminoacid_length(self):
        return len(self.sequence)


    def add_ribosome(self, ribosome):
        """
        By definition this reaction will be catalyzed by a ribosome
        :param ribosome:
        :type ribosome: pytfa.me.Ribosome
        :return:
        """
        self.enzyme = ribosome


class ProteinComplexation(Reaction):
    """
    Describes the assembly of peptides into an enzyme
    """
    def __init__(self, *args, **kwargs):
        Reaction.__init__(self, *args, **kwargs)
        self.enzyme = None