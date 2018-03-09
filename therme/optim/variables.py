# -*- coding: utf-8 -*-
"""
.. module:: thermome
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Variables declarations

"""

from pytfa.optim.variables import GenericVariable, BinaryVariable, \
    ReactionVariable, get_binary_type

class ModelVariable(GenericVariable):
    """
    Class to represent a variable attached to the model
    """

    def __init__(self, model, id_, **kwargs):
        if not 'lb' in kwargs:
            kwargs['lb'] = 0
        GenericVariable.__init__(self,
                                 id_= id_,
                                 model=model,
                                 **kwargs)


class GrowthRate(ModelVariable):
    """
    Class to represent a growth rate
    """

    def __init__(self, model, **kwargs):
        id_ = 'mu'
        ModelVariable.__init__(self, model=model, id_= id_, **kwargs)

    prefix = 'MU_'


class GrowthActivation(ModelVariable, BinaryVariable):
    """
    Class to represent a binary growth rate range activation in ME2 MILP
    """

    def __init__(self, model, id_, **kwargs):
        self.ix = int(id_)
        BinaryVariable.__init__(self, id_= id_, model=model, **kwargs)

    prefix = 'GA_'


class EnzymeVariable(GenericVariable):
    """
    Class to represent a enzyme variable
    """

    prefix = 'EZ_'

    def __init__(self, enzyme, **kwargs):
        self.enzyme = enzyme
        model = enzyme.model


        if not 'lb' in kwargs:
            kwargs['lb'] = 0
        if not 'ub' in kwargs:
            kwargs['ub'] = model.max_enzyme_concentration


        GenericVariable.__init__(self, id_=self.id, model=model,
                                 **kwargs)

    @property
    def id(self):
        return self.enzyme.id

    @property
    def model(self):
        return self.enzyme.model

class GeneVariable(GenericVariable):
    """
    Class to represent a gene variable
    """

    prefix = 'GV_'

    def __init__(self, gene, **kwargs):
        self.gene = gene
        model = gene.model


        if not 'lb' in kwargs:
            kwargs['lb'] = 0
        if not 'ub' in kwargs:
            kwargs['ub'] = model.max_enzyme_concentration

        GenericVariable.__init__(self, id_=self.id, model=model,
                                 **kwargs)

    @property
    def id(self):
        return self.gene.id

    @property
    def model(self):
        return self.gene.model



class mRNAVariable(GeneVariable):
    """
    Class to represent a mRNA concentration
    """

    prefix = 'MR_'

class ForwardEnzyme(EnzymeVariable):
    """
    Represents assignment of an enzyme the a forward reaction flux
    """
    prefix = 'FE_'

class BackwardEnzyme(EnzymeVariable):
    """
    Represents assignment of an enzyme the a backward reaction flux
    """
    prefix = 'BE_'


class LinearizationVariable(ModelVariable):
    """
    Class to represent the product mu*[E] when performin linearization of the
    model
    """
    prefix = 'LZ_'


class RibosomeUsage(GeneVariable):
    """
    Class to represent the ribosomes that are assigned to producing the enzyme
    for a reaction
    """
    prefix = 'RP_'


class RNAPUsage(GeneVariable):
    """
    Class to represent the ribosomes that are assigned to producing the enzyme
    for a reaction
    """
    prefix = 'RM_'


class FreeRibosomes(EnzymeVariable):
    """
    Class to represent the ribosomes that are affected to producing the enzyme
    for a reaction
    """
    prefix = 'RF_'


class CatalyticActivator(ReactionVariable,BinaryVariable):
    """
    Class to represent a binary variable that activates a catalytic constraint
    or relaxes it
    """
    def __init__(self, reaction, **kwargs):
        ReactionVariable.__init__(self, reaction=reaction,
                                  type = get_binary_type(),
                                  **kwargs)

    prefix = 'CA_'

class BinaryActivator(ModelVariable, BinaryVariable):
    """
    Class to represent a binary variable that activates with growth levels
    """
    def __init__(self, model, id_, **kwargs):
        ModelVariable.__init__(self,
                               model,
                               id_,
                               type = get_binary_type(),
                               **kwargs)

    prefix = 'LA_'

class InterpolationVariable(ModelVariable):
    """
    Represents a variable that is interpolated
    """
    prefix = 'IV_'