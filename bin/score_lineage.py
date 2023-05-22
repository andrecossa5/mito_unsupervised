#!/usr/bin/python

# Score LINEAGE output

import os
import sys
import pandas as pd
import pickle
from sklearn.metrics import normalized_mutual_info_score
from mito_utils.utils import *
from mito_utils.clustering import *
from mito_utils.preprocessing import *

# Path
sample = sys.argv[1]
path_ground_truth = sys.argv[2]
path_labels = sys.argv[3]
min_cell_number = sys.argv[3]

# Read
gt = pd.read_csv(path_ground_truth, index_col=0)['GBC']
labels = pd.read_csv(path_labels, index_col=0)['label']

# Score
L = []
ari = custom_ARI(labels, gt)
nmi = normalized_mutual_info_score(labels, gt)
     
L.append({
    'model' : 'LINEAGE',
    'ARI' : ari,
    'NMI' : nmi,
    'sample' : sample,
    'filtering': 'custom_lineage',
    'dimred' : 'no_dimred',
    'min_cell_number' : min_cell_number,
    'min_cov_treshold' : 50,
    'ncells_sample' : gt.size,
    'n_clones_analyzed' : gt.unique().size,
    'n_clones_inferred' : np.unique(labels).size,
    'n_variants' : 'custom_lineage'
})
    
# Write
df = pd.DataFrame(L).sort_values('NMI')
d = {'df_performance' : df, 'labels' : labels}
with open(f'out_LINEAGE_{sample}_{min_cell_number}.pickle', 'wb') as f:
    pickle.dump(d, f)