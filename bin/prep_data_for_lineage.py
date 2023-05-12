#!/usr/bin/python

# Prep AFM for LINEAGE

import os
import sys
import scanpy as sc
from mito_utils.utils import *
from mito_utils.preprocessing import *

# Path
path_data = sys.argv[1]
sample = sys.argv[2]

# Read
afm = read_one_sample(path_data, sample=sample)

# Write
pd.DataFrame(
    afm.X, 
    index=afm.obs_names,
    columns=afm.var_names,
).to_csv('formatted.csv')

# meta
afm.obs.to_csv('meta.csv')