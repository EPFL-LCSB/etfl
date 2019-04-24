# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor:: ETFL team

Utilities to create a small model from a 1 reaction model in FBA

"""
import pytest

from etfl.tests.small_model import create_etfl_model
from etfl.io.json import json_dumps_model, json_loads_model, save_json_model, load_json_model
from etfl.io.dict import model_from_dict, model_to_dict
# Dependency help from
# https://stackoverflow.com/questions/49238725/chaining-tests-and-passing-an-object-from-one-test-to-another

# Caching help from
# https://stackoverflow.com/questions/22768976/how-to-share-a-variable-across-modules-for-all-tests-in-py-test/50611422#50611422


@pytest.fixture(autouse=True)
def init_cache(request):
    m = request.config.cache.get('small_etfl', None)
    m = create_etfl_model(has_thermo=False,
                          has_neidhardt=False,
                          n_mu_bins = 4,
                          optimize = False)
    request.config.cache.set('small_etfl', json_dumps_model(m))

@pytest.mark.dependency(name='test_etfl')
def test_etfl(request):
    model = json_loads_model(request.config.cache.get('small_etfl',None))
    assert model is not None
    model.optimize()

savepath = './test_model.json'

@pytest.mark.dependency()
@pytest.mark.dependency(depends=['test_etfl'])
def test_write(request):
    model = json_loads_model(request.config.cache.get('small_etfl',None))
    assert model is not None
    save_json_model(model, savepath)


@pytest.mark.dependency(depends=['test_write'])
def test_read():
    m = load_json_model(savepath)
    m.optimize()


# Optim configs
from etfl.optim.config import standard_solver_config, gene_ko_config, \
    growth_uptake_config

@pytest.mark.dependency(depends=['test_etfl'])
def test_configs(request):
    model = json_loads_model(request.config.cache.get('small_etfl', None))
    standard_solver_config(model)
    gene_ko_config(model)
    growth_uptake_config(model)