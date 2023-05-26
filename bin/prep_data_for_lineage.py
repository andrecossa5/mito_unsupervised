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
min_cell_number = sys.argv[3]

# Read
afm = read_one_sample(path_data, sample=sample)

# Filter good quality cells
a_cells, _ = filter_cells_and_vars(
    afm,
    sample=sample,
    filtering='LINEAGE_prep',
    min_cell_number=int(min_cell_number),
    min_cov_treshold=50
)

# Write AFM and meta
pd.DataFrame(
    a_cells.X, 
    index=a_cells.obs_names,
    columns=a_cells.var_names,
).to_csv('formatted.csv')
a_cells.obs.to_csv('meta.csv')


##