#!/usr/bin/python

# Prep json for CloneTracer

import os
import json
import sys
import scanpy as sc
from mito_utils.utils import *
from mito_utils.preprocessing import *

# Path
# path_data = sys.argv[1]
# sample = sys.argv[2]
# min_cell_number = sys.argv[3]

path_data = '/Users/IEO5505/Desktop/mito_bench/data/'
sample = 'MDA_PT'
min_cell_number = 10

# Read
afm = read_one_sample(path_data, sample=sample, with_GBC=True)

# Filter good quality cells
_, a = filter_cells_and_vars(
    afm,
    sample=sample,
    filtering='MQuad',
    min_cell_number=int(min_cell_number),
    min_cov_treshold=50,
    path_=os.getcwd()
)

# Extract filtered feature matrix, format and reduce with UMAP
a = nans_as_zeros(a) # For sklearn APIs compatibility

# Get parallel matrices
AD, DP, _ = get_AD_DP(a, to='csc')

# Prep dict
d = {
    'M' : AD.A.T.tolist(),
    'N' : DP.A.T.tolist(),
    'mut_names' : a.var_names.to_list(),
    'cell_barcode' : {'V1':a.obs_names.to_list()},
    'mut_type' : [ 2 for _ in a.var_names ]
}

# Write as json
with open(os.path.join(path_data, 'MDA_PT_prova.json'), "w") as f:
    json.dump(d, f)
    

##