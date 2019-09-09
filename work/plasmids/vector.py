# -*- coding: utf-8 -*-
"""
.. module:: ETFL
   :platform: Unix, Windows
   :synopsis: Models vectors in ETFL models

.. moduleauthor:: ETFL team

Model RNAP limitation in case of vector addition to a host
"""

from etfl.io.json import load_json_model
from etfl.core.vector import TransModel, make_plasmid
from etfl.core.genes import ExpressedGene

from utils import read_seq

# E. coli

model = load_json_model('models/SlackModel '
                    'iJO1366_vETFL_tp_441_enz_128_bins__20190620_144322.json')

# Plasmid pZS*-13S-ald-adh
# US Patent 20,140,371,417.
# Used in Andreozzi, Stefano, et al.
# "Identification of metabolic engineering targets for the enhancement of
# 1,4-butanediol production in recombinant E. coli using large-scale
# kinetic models."
# Metabolic engineering 35 (2016): 148-159.

# (too complicated)

# 2-gene BDO plasmid
# Reshamwala, Shamlan MS, Shalini S. Deb, and Arvind M. Lali.
# "A shortened, two-enzyme pathway for 2, 3-butanediol production in
# Escherichia coli."
# Journal of industrial microbiology & biotechnology 44.9 (2017): 1273-1277.

ALS = ExpressedGene(id = 'eba:als',
                    name = 'acetolactate synthase',
                    sequence = read_seq('data/ALS.txt'))
AR = ExpressedGene (id = 'eba:ar',
                    name = 'acetoin reductase',
                    sequence = read_seq('data/AR.txt'))
fwd_als= 'GGAATTCCATATGAACAGTGAGAAACAGTC'
rev_als= 'CCGAGCTCTTACAAAATCTGGCTGAGAT'
fwd_ar = 'GGAATTCCATATGCAAAAAGTTGCTCTCGTAAC'
rev_ar = 'CCGAGCTCTTAGTTGAACACCATCCCAC'

pET = read_seq('data/pET.txt')

gene_list = [ALS,AR]

plasmid_seq = pET + fwd_als + ALS.seq + rev_als + fwd_ar + AR.seq + rev_ar

vector = make_plasmid()

