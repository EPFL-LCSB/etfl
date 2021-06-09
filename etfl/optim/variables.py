# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

Variables declarations

"""

from pytfa.optim.variables import GenericVariable, BinaryVariable, \
    ReactionVariable, ModelVariable, GeneVariable, get_binary_type




class GrowthRate(ModelVariable):
    """
    Class to represent a growth rate
    """

    def __init__(self, model, **kwargs):
        ModelVariable.__init__(self, model=model, **kwargs)

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
        model = enzyme.model


        if not 'lb' in kwargs:
            kwargs['lb'] = 0
        if not 'ub' in kwargs:
            kwargs['ub'] = model.max_enzyme_concentration


        GenericVariable.__init__(self, model=model, hook=enzyme,
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



class mRNAVariable(GeneVariable):
    """
    Class to represent a mRNA concentration
    """

    prefix = 'MR_'

class rRNAVariable(GeneVariable):
    """
    Class to represent a mRNA concentration
    """

    prefix = 'RR_'

class tRNAVariable(ModelVariable):
    """
    Class to represent a tRNA concentration
    """

    prefix = 'TR_'

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


class DNAVariable(ModelVariable):
    """
    Class to represent DNA in the model
    """
    prefix = 'DN_'


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


class FreeEnzyme(EnzymeVariable):
    """
    Class to represent the ribosomes that are affected to producing the enzyme
    for a reaction
    """
    prefix = 'EF_'


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

class EnzymeRef(EnzymeVariable):
    """
    Represents a reference enzyme concentration - for example in dETFL
    """

    prefix = 'EZ0_'

class mRNARef(mRNAVariable):
    """
    Represents a reference enzyme concentration - for example in dETFL
    """

    prefix = 'MR0_'
    
class LipidVariable(ModelVariable):
    """
    Class to represent lipid in the model
    """
    prefix = 'LIP_'
    
class CarbohydrateVariable(ModelVariable):
    """
    Class to represent carbohydrate in the model
    """
    prefix = 'CAR_'
    
class IonVariable(ModelVariable):
    """
    Class to represent ion in the model
    """
    prefix = 'ION_'