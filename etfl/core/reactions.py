# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

ME-related Reaction subclasses and methods definition


"""
from cobra import Reaction, DictList
from _collections import defaultdict


class ExpressionReaction(Reaction):
    def __init__(self, scaled, **kwargs):
        Reaction.__init__(self, **kwargs)
        self._scaled = scaled

    @classmethod
    def from_reaction(cls,reaction, scaled = False, **kwargs):
        """
        This method clones a cobra.Reaction object into a expression-related
        type of reaction

        :param reaction: the reaction to reproduce
        :return: an EnzymaticReaction object
        """
        new =  cls( id = reaction.id,
                    name= reaction.name,
                    subsystem= reaction.subsystem,
                    lower_bound= reaction.lower_bound,
                    upper_bound= reaction.upper_bound,
                    scaled = scaled,
                    **kwargs)

        # new._model = reaction._model
        # new.notes = reaction.notes
        # We need not rescale the initial metabolites
        new.add_metabolites(reaction.metabolites, rescale=False)
        new.gene_reaction_rule = reaction.gene_reaction_rule
        return new

    def add_metabolites(self, metabolites, rescale = True, **kwargs):
        """
        We need to override this method if the reaction is scaled

        v_hat = v/vmax

        dM/dt = n1*v1 + ...

        dM/dt = n1*vmax1 * v1_hat + ...

        :param metabolites:
        :return:
        """

        if not hasattr(self, '_scaled') or not rescale or not self._scaled:
            Reaction.add_metabolites(self, metabolites, **kwargs)
        else:
            Reaction.add_metabolites(self, {k:v*self.scaling_factor
                                        for k,v in metabolites.items()},
                                     **kwargs)

    @property
    def scaling_factor(self):
        return 1

    @property
    def net(self):
        return self.scaling_factor * self.scaled_net

    @property
    def scaled_net(self):
        return self.forward_variable - self.reverse_variable

class EnzymaticReaction(ExpressionReaction):
    """
    Subclass to describe reactions that are catalyzed by an enzyme.
    """
    def __init__(self, enzymes = None, scaled = False, *args, **kwargs):
        ExpressionReaction.__init__(self, scaled, *args, **kwargs)
        self.enzymes = DictList()
        if enzymes:
            self.add_enzymes(enzymes)

    def add_enzymes(self, enzymes):
        """`
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

    @property
    def scaling_factor(self):
        return self.enzyme.kcat_max * self.enzyme.scaling_factor
        # return 1


class TranscriptionReaction(EnzymaticReaction):
    """
    Class describing transcription - Assembly of amino acids into peptides
    """

    def __init__(self, id, name, gene_id, enzymes, **kwargs):
        EnzymaticReaction.__init__(self,
                                   id=id,
                                   name=name,
                                   enzymes=enzymes,
                                   **kwargs)
        self._gene_id = gene_id

    # The number of amino acids is in the stoichiometry

    @property
    def gene(self):
        return self.model.genes.get_by_id(self._gene_id)

    @property
    def nucleotide_length(self):
        return len(self.gene.rna)


    def add_rnap(self, rnap):
        """
        By definition this reaction will be catalyzed by RNA polymerase
        :param ribosome:
        :type ribosome: pytfa.me.RNAPolymerase
        :return:
        """
        self.enzymes = rnap

    @property
    def scaling_factor(self):
        return self.enzymes[0].kcat_fwd * self.enzymes[0].scaling_factor \
                / self.nucleotide_length
        # return 1

class TranslationReaction(EnzymaticReaction):
    """
    Class describing translation - Assembly of amino acids into peptides
    """

    def __init__(self, id, name, gene_id, enzymes, trna_stoich=None, **kwargs):
        EnzymaticReaction.__init__(self,
                                   id=id,
                                   name=name,
                                   enzymes=enzymes,
                                   **kwargs)
        self._gene_id = gene_id
        if trna_stoich is None:
            self.trna_stoich = defaultdict(int)
        else:
            self.trna_stoich = trna_stoich

    @property
    def gene(self):
        return self.model.genes.get_by_id(self._gene_id)

    # The number of amino acids is in the stoichiometry
    @property
    def aminoacid_length(self):
        return len(self.gene.peptide)

    def add_peptide(self, peptide):
        """
        According to the scaling rules, the coefficient of the scaled translation
        reaction for the peptide balance is 1:

        dPep/dt = v_tsl     - sum(ηj * vj_asm) = 0
                  v_tsl_hat - sum(ηj * L_aa/(krib * R_max) * kdegj * Ej_max * vj_asm_max)

        :param peptide:
        :return:
        """

        self.add_metabolites({peptide:1}, rescale=True)


    def add_ribosome(self, ribosome):
        """
        By definition this reaction will be catalyzed by a ribosome
        :param ribosome:
        :type ribosome: pytfa.me.Ribosome
        :return:
        """
        self.enzymes = ribosome



    @property
    def scaling_factor(self):
        return self.enzymes[0].kcat_fwd * self.enzymes[0].scaling_factor \
                / self.aminoacid_length
        # return 1


class ProteinComplexation(ExpressionReaction):
    """
    Describes the assembly of peptides into an enzyme
    """
    def __init__(self, target, *args, **kwargs):
        ExpressionReaction.__init__(self, *args, **kwargs)
        self.enzymes = None
        self.target = target
        self.target.complexation = self

    @property
    def scaling_factor(self):
        # return self.model.mu_max * self.target.scaling_factor
        # return self.target.kdeg * self.target.scaling_factor
        return 1

    def add_peptides(self, peptides):
        """
        /!\ Reaction must belong to a model

        According to the scaling rules, the coefficient of the scaled complexation
        reaction for the peptide balance is L_aa/(krib * R_max):

        dPep/dt = v_tsl     - sum(ηj * vj_asm) = 0
                  v_tsl_hat - sum(ηj * L_aa/(krib * R_max) * kdegj * Ej_max * vj_asm_max)


        :param peptides: dict(Peptide: int)
        :return:
        """

        self.add_metabolites({p:s for p,s in peptides.items()}, rescale=True)


class DegradationReaction(ExpressionReaction):
    """
    Describes the degradation of macromolecules
    """
    def __init__(self, macromolecule, *args, **kwargs):
        ExpressionReaction.__init__(self, *args, **kwargs)
        self.enzymes = None
        self.macromolecule = macromolecule
        self.macromolecule.degradation = self


    @property
    def scaling_factor(self):
        # return self.model.mu_max * self.macromolecule.scaling_factor
        return self.macromolecule.kdeg * self.macromolecule.scaling_factor
        # return 1

class DNAFormation(ExpressionReaction):
    """
    Describes the assembly of NTPs into DNA
    """
    def __init__(self, dna, mu_sigma=1, *args, **kwargs):
        ExpressionReaction.__init__(self, *args, **kwargs)
        self.dna = dna
        # mu_sigma is a scaling factor ~ mu_max (same homogeneity)
        self.mu_sigma = mu_sigma

    @property
    def scaling_factor(self):
        return self.mu_sigma * self.dna.scaling_factor