# -*- coding: utf-8 -*-
"""
.. module:: etfl
   :platform: Unix, Windows
   :synopsis: Expression and thermodynamics-based models

.. moduleauthor:: Pierre Salvy

JSON serialization
"""

import json

from .dict import model_from_dict, model_to_dict
from pytfa.io.json import MyEncoder, check_json_extension

def save_json_model(model, filepath):
    """
    Saves the model as a JSON file

    :param model:
    :param filepath:
    :return:
    """

    filepath = check_json_extension(filepath)
    obj = model_to_dict(model)

    with open(filepath, 'w') as fid:
        json.dump(obj, fid, cls=MyEncoder)


def load_json_model(filepath, solver = None):
    """
    Loads a model from a JSON file

    :param filepath:
    :param solver:
    :return:
    """

    filepath = check_json_extension(filepath)
    with open(filepath, 'r') as fid:
        obj = json.load(fid)

    model = model_from_dict(obj, solver=solver)
    return model

def json_dumps_model(model):
    """
    Returns a JSON dump as a string

    :param model:
    :return:
    """

    obj = model_to_dict(model)

    return json.dumps(obj,cls=MyEncoder)


def json_loads_model(s):
    """
    Loads a model from a string JSON dump

    :param s: JSON string
    :return:
    """
    obj = json.loads(s)
    return model_from_dict(obj)