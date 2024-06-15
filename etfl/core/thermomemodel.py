"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

Fusion for Thermo and Me Models

"""

from pytfa.core import LCSBModel
from .memodel import MEModel
from pytfa.thermo import ThermoModel, std, get_debye_huckel_b

from cobra import Model, DictList

from pytfa.utils.logger import get_bistream_logger
from pytfa.utils import numerics
from pytfa.optim.utils import copy_solver_configuration


import optlang

BIGM = numerics.BIGM
BIGM_THERMO = numerics.BIGM_THERMO
BIGM_DG = numerics.BIGM_DG
BIGM_P = numerics.BIGM_P
EPSILON = numerics.EPSILON
MAX_STOICH = 10

class ThermoMEModel(MEModel, ThermoModel):

    def __init__(self,
                 thermo_data,
                 model=Model(),
                 name = None,
                 growth_reaction='',
                 mu = None, mu_error = 0,
                 mu_range = None, n_mu_bins = 1,
                 big_M = 1000,
                 temperature=std.TEMPERATURE_0,
                 min_ph=std.MIN_PH,
                 max_ph=std.MAX_PH,
                 prot_scaling = 1000,
                 mrna_scaling = None, ):

        if name is None:
            name = 'ETFL_' + name if name else 'ETFL model'

        LCSBModel.__init__(self, model, name)

        ###############
        #   ME part   #
        ###############

        self.logger = get_bistream_logger('ME model' + str(self.name))
        self.parent = model
        # if model is not None:
        #     self.sanitize_varnames()

        self.init_etfl(big_M, growth_reaction, mu_range, n_mu_bins, name)

        ###############
        # Thermo part #
        ###############

        self.TEMPERATURE = temperature
        self.thermo_data = thermo_data
        self.thermo_unit = thermo_data['units']
        self.reaction_cues_data = thermo_data['cues']
        self.compounds_data = thermo_data['metabolites']
        self.Debye_Huckel_B = get_debye_huckel_b(temperature)
        self.parent = model

        self.logger = get_bistream_logger('thermomodel_' + str(self.name))

        # Compute internal values to adapt the the thermo_unit provided
        if self.thermo_unit == "kJ/mol":
            self.GAS_CONSTANT = 8.314472 / 1000  # kJ/(K mol)
            self.Adjustment = 1
        else:
            self.GAS_CONSTANT = 1.9858775 / 1000  # Kcal/(K mol)
            self.Adjustment = 4.184

        self.RT = self.GAS_CONSTANT * self.TEMPERATURE

        # CONSTANTS
        self.MAX_pH = max_ph
        self.MIN_pH = min_ph

        self.logger.info('# Model initialized with units {} and temperature {} K'  \
                    .format(self.thermo_unit, self.TEMPERATURE))


    def print_info(self):
        LCSBModel.print_info(self)
        ThermoModel.print_info(self, specific=True)
        MEModel.print_info(self, specific=True)

    def __deepcopy__(self, memodict={}):
        return self.copy()

    def copy(self):

        from ..io.dict import model_from_dict, model_to_dict
        dictmodel = model_to_dict(self)
        new = model_from_dict(dictmodel)

        copy_solver_configuration(self, new)

        return new