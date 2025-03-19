import os
import sys
import itertools
import pandas as pd
import numpy as np

script_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'data', 'aa_physiochem_metrics', 'rough'))
sys.path.append(script_dir)

# Get AA distance/similarity matrices for Miyata, EMPAR, CSW, and EX.

# Note that I need to bump to >= Python 3.10 to import this module.
# This Python file was downloaded from commit 5c428e558f0c61a65549ac3fc62a6a8024a06fc0
# of the https://github.com/MICS-Lab/pard/ GitHub repository.
import _raw_python_dictionaries

amino_acids = ["F", "L", "S", "Y", "C", "W", "P", "H", "Q", "R", "I", "M",
               "T", "N", "K", "V", "A", "D", "E", "G"]

pairs = list(itertools.combinations(amino_acids, 2))

empty_df = pd.DataFrame(np.zeros((20, 20)), index=amino_acids, columns=amino_acids)

# Miyata (This is the only distance - others are similarities)
miyata_dist = empty_df.copy(deep=True)
MIYATA_DICT = _raw_python_dictionaries.MIYATA_DICT
for pair in pairs:
    miyata_dist.loc[pair[0], pair[1]] = MIYATA_DICT[(pair[0], pair[1])]
    miyata_dist.loc[pair[1], pair[0]] = MIYATA_DICT[(pair[1], pair[0])]
miyata_dist.to_csv('/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/distances/miyata_orig.tsv',
                   index=True, header=True, sep='\t')

# EX
ex_sim = empty_df.copy(deep=True)
EX_DICT = _raw_python_dictionaries.SYMMETRIC_EXPERIMENTAL_EXCHANGEABILITY_DICT
for pair in pairs:
    ex_sim.loc[pair[0], pair[1]] = EX_DICT[(pair[0], pair[1])]
    ex_sim.loc[pair[1], pair[0]] = EX_DICT[(pair[1], pair[0])]

# Set diagonal to NaN.
np.fill_diagonal(ex_sim.values, np.nan)
ex_sim.to_csv('/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/similarities/ex.tsv',
               index=True, header=True, sep='\t')

# EMPAR
empar_sim = empty_df.copy(deep=True)
EMPAR_DICT = _raw_python_dictionaries.EMPAR_DICT
for pair in pairs:
    empar_sim.loc[pair[0], pair[1]] = EMPAR_DICT[(pair[0], pair[1])]
    empar_sim.loc[pair[1], pair[0]] = EMPAR_DICT[(pair[1], pair[0])]

# Set diagonal to 16.
np.fill_diagonal(empar_sim.values, 16)
empar_sim.to_csv('/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/similarities/empar.tsv',
                  index=True, header=True, sep='\t')

# CSW
csw_sim = empty_df.copy(deep=True)
CSW_DICT = _raw_python_dictionaries.KOLASKAR_DICT
for pair in pairs:
    csw_sim.loc[pair[0], pair[1]] = CSW_DICT[(pair[0], pair[1])]
    csw_sim.loc[pair[1], pair[0]] = CSW_DICT[(pair[1], pair[0])]
# Set diagonal to 10.
np.fill_diagonal(csw_sim.values, 10)
csw_sim.to_csv('/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/similarities/csw.tsv',
                index=True, header=True, sep='\t')
