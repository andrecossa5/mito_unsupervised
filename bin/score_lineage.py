#!/usr/bin/python

# Score LINEAGE output

import os
import sys
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score
from mito_utils.utils import *
from mito_utils.clustering import *
from mito_utils.preprocessing import *

# Path
sample = sys.argv[1]
path_ground_truth = sys.argv[2]
path_labels = sys.argv[3]

# Read
gt = pd.read_csv(path_ground_truth, index_col=0)
labels = pd.read_csv(path_labels, index_col=0)

# Write
L = []
ari = custom_ARI(labels['label'], gt['GBC'])
nmi = normalized_mutual_info_score(labels['label'], gt['GBC'])
     
L.append({
    'ARI' : ari,
    'NMI' : nmi,
    'sample' : sample
})

# Write 
(
    pd.DataFrame(L)
    .sort_values('ARI')
    .to_csv(f'out_{sample}.csv')
)