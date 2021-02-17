# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

Constraints declarations

"""
from pytfa.optim import ReactionConstraint, GenericConstraint, \
    ModelConstraint, GeneConstraint


class CatalyticConstraint(ReactionConstraint):
    """
    Class to represent a enzymatic constraint
    """

    prefix = 'CC_'


class ForwardCatalyticConstraint(ReactionConstraint):
    """
    Class to represent a enzymatic constraint
    """

    prefix = 'FC_'


class BackwardCatalyticConstraint(ReactionConstraint):
    """
    Class to represent a enzymatic constraint
    """

    prefix = 'BC_'


class EnzymeConstraint(GenericConstraint):
    """
    Class to represent a variable attached to a enzyme
    """

    def __init__(self, enzyme, expr, **kwargs):
        model = enzyme.model

        GenericConstraint.__init__(self,
                                   expr=expr,
                                   model=model,
                                   hook=enzyme,
                                   **kwargs)

    @property
    def enzyme(self):
        return self.hook

    @property
    def id(self):
        return self.enzyme.id

    @property
    def model(self):
        return self.enzyme.model

    prefix = 'EZ_'


class EnzymeMassBalance(EnzymeConstraint):
    """
    Class to represent a enzymatic mass balance constraint
    """

    prefix = 'EB_'

class mRNAMassBalance(GeneConstraint):
    """
    Class to represent a mRNA mass balance constraint
    """

    prefix = 'MB_'

class rRNAMassBalance(GeneConstraint):
    """
    Class to represent a mRNA mass balance constraint
    """

    prefix = 'RB_'

class tRNAMassBalance(ModelConstraint):
    """
    Class to represent a tRNA mass balance constraint
    """

    prefix = 'TB_'

class DNAMassBalance(ModelConstraint):
    """
    Class to represent a DNA mass balance constraint
    """

    prefix = 'DB_'


class SynthesisConstraint(ReactionConstraint):
    """
    Class to represent a Translation constraint
    """

    prefix = 'TR_'


class GrowthCoupling(ReactionConstraint):
    """
    Class to represent a growth capacity constraint
    """

    prefix = 'GC_'


class TotalCapacity(ModelConstraint):
    """
    Class to represent the total capacity of constraint of a species, e.g
    Ribosome or RNA
    """

    prefix = 'TC_'


class TotalEnzyme(TotalCapacity):
    """
    Class to represent the total amount of an enzyme species, forwards and backwards
    """

    prefix = 'TE_'


class ExpressionCoupling(GeneConstraint):
    """
    Add the coupling between mRNA availability and ribosome charging
    The number of ribosomes assigned to a mRNA species is lower than
    the number of such mRNA times the max number of ribosomes that can sit
    on the mRNA:
    [RPi] <= loadmax_i*[mRNAi]
    """

    prefix = 'EX_'
    
class MinimalCoupling(GeneConstraint):
    """
    Add the minimal activity of ribosome based on the availability of mRNA.
    We modeled it as a fraction of the maximum loadmax and the fraction depends
    on the affinity of ribosome to the mRNA:
    [RPi] >= Fraction*loadmax_i*[mRNAi]
    """

    prefix = 'MC_'


class RNAPAllocation(GeneConstraint):
    """
    Add the coupling between DNA availability and RNAP charging
    The number of RNAP assigned to a gene locus is lower than
    the number of such loci times the max number of RNAP that can sit
    on the locus:
    [RNAPi] <= loadmax_i*[# of loci]*[DNA]
    """

    prefix = 'RA_'
    
class MinimalAllocation(GeneConstraint):
    """
    Add the minimal activity of RNAP based on the availability of gene.
    We modeled it as a fraction of the maximum loadmax and the fraction depends
    on the affinity of RNAP to the gene, i.e. the strength of the promoter:
    [RPi] >= Fraction*loadmax_i*[mRNAi]
    """

    prefix = 'MA_'


class EnzymeRatio(EnzymeConstraint):
    """
    Represents the availability of free enzymes, e.g ribosomes (non bound)
    R_free = 0.2*R_total
    """

    prefix = 'ER_'

class RibosomeRatio(EnzymeRatio):
    """
    (Legacy) represents the availability of free ribosomes, e.g ribosomes (non bound)
    R_free = 0.2*R_total
    """

    prefix = 'ER_'

class EnzymeDegradation(EnzymeConstraint):
    """
    v_deg = k_deg [E]
    """

    prefix = "ED_"

class mRNADegradation(GeneConstraint):
    """
    v_deg = k_deg [mRNA]
    """

    prefix = "MD_"


class GrowthChoice(ModelConstraint):
    """
    Class to represent a variable attached to a reaction
    """

    prefix = 'GR_'


class LinearizationConstraint(ModelConstraint):
    """
    Class to represent a variable attached to a reaction
    """
    @staticmethod
    def from_constraints(cons, model):
        return LinearizationConstraint(
            name = cons.name,
            expr = cons.expr,
            model = model,
            ub = cons.ub,
            lb = cons.lb,
        )

    prefix = 'LC_'

class SOS1Constraint(ModelConstraint):
    """
    Class to represent SOS 1 constraint
    """

    prefix = 'S1_'

class InterpolationConstraint(ModelConstraint):
    """
    Class to represent an interpolation constraint
    """

    prefix = "IC_"


class EnzymeDeltaPos(EnzymeConstraint):
    """
    Represents a positive enzyme concentration variation for dETFL
    """

    prefix = 'dEP_'

class EnzymeDeltaNeg(EnzymeConstraint):
    """
    Represents a negative enzyme concentration variation for dETFL
    """

    prefix = 'dEN_'


class mRNADeltaPos(GeneConstraint):
    """
    Represents a positive mRNA concentration variation for dETFL
    """

    prefix = 'dMP_'

class mRNADeltaNeg(GeneConstraint):
    """
    Represents a negative mRNA concentration variation for dETFL
    """

    prefix = 'dMN_'
    
class ConstantAllocation(ModelConstraint):
    """
    Represents a similar share to FBA for RNA and protein
    """
    
    prefix = 'CL_'
    
class LipidMassBalance(ModelConstraint):
    """
    Class to represent a lipid mass balance constraint
    """

    prefix = 'LB_'
    
class CarbohydrateMassBalance(ModelConstraint):
    """
    Class to represent a carbohydrate mass balance constraint
    """

    prefix = 'CB_'
    
class IonMassBalance(ModelConstraint):
    """
    Class to represent a ion mass balance constraint
    """

    prefix = 'IB_'