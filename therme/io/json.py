# -*- coding: utf-8 -*-
"""
.. module:: therme
   :platform: Unix, Windows
   :synopsis: Expression and thermodynamics-based models

.. moduleauthor:: Pierre Salvy

JSON serialization
"""

import json

from .dict import model_from_dict, model_to_dict
from pytfa.io.json import MyEncoder, check_json_extension

def save_json_model(model, filepath):

    filepath = check_json_extension(filepath)
    obj = model_to_dict(model)

    with open(filepath, 'w') as fid:
        json.dump(obj, fid, cls=MyEncoder)


def load_json_model(filepath):

    filepath = check_json_extension(filepath)
    with open(filepath, 'r') as fid:
        obj = json.load(fid)

    model = model_from_dict(obj)
    return model
