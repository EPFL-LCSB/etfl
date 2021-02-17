# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 18:56:43 2020

@author: Omid
"""
### carbohydrate.py
from .macromolecule import Macromolecule
from ..optim.variables import CarbohydrateVariable

class Carbohydrate(Macromolecule):
    def __init__(self, composition, mass_ratio, id='Carbohydrate', kdeg=0, *args, **kwargs):
        Macromolecule.__init__(self, id = id, kdeg=kdeg, *args, **kwargs)
        self.composition = composition
        self._MW = mass_ratio

    def init_variable(self, queue=False):
        """
        Attach a carbohydrateVariable object to the Species. Needs to have the object
        attached to a model

        :return:
        """
        self._internal_variable = self.model.add_variable(kind=CarbohydrateVariable,
                                    hook=self.model,
                                    id_=self.id,
                                    lb=0,
                                    ub=1)
    @property
    def molecular_weight(self):
        """
        To keep consistency

        :return:
        """
        return self._MW # the mass ratio in FBA biomass