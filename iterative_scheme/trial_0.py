"""
Prep MDA_PT data for iterative_scheme.
"""

# Code
import os
import pickle
from kneed import KneeLocator
from vireoSNP import BinomMixtureVB
from sklearn.metrics import normalized_mutual_info_score

from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.clustering import *
from mito_utils.iterative import *


##


# Set paths
path_main = '/Users/IEO5505/Desktop/mito_bench/'
path_data = os.path.join(path_main, 'PT_subsampled') # Here subsampled
path_output = os.path.join(path_main, 'results', 'unsupervised_clones', 'output')
path_viz = os.path.join(path_main, 'results', 'unsupervised_clones', 'visualization')
path_tmp = os.path.join(path_main, 'results', 'unsupervised_clones', 'downstream_files')

##

# AFM first step, no filtering of low numebr of cells clones + variants filtering
afm = read_one_sample(path_data, 'PT_subsampled', with_GBC=True)
a_cells, a = filter_cells_and_vars(afm, sample='PT_subsampled', filtering='MQuad',
                                   min_cov_treshold=50, nproc=8, path_=os.getcwd()) # Do not exclude rare clones here!

##

# Read and subset cells and vars

# Here we go
i = 0
go = True
with_ex = {}
variants = {}

##

# Read saved info for iteration 0
cells = pd.read_csv(os.path.join(path_tmp, 'MDA_PT_filtered_cells'), index_col=0).values.flatten().tolist()
variants = pd.read_csv(os.path.join(path_tmp, 'MDA_PT_original_variants'), index_col=0).values.flatten().tolist()
a_cells, a = filter_cells_and_vars(afm, sample='MDA_PT', cells=cells, variants=variants)

# Iteration 0:
a.obs[f'it_{i}'] = vireo_wrapper(a)

a


afm.obs[f'it_{i}'].value_counts().to_dict()
exclusive_variants = { 
    c : rank_clone_variants(afm, var=f'it_{i}', group=c, 
                            mode = 'perc_ratio',
                            log2_min_perc_ratio=4).index.to_list() \
    for c in afm.obs[f'it_{i}'].unique() 
}
wi, wo = test_partitions_variants(exclusive_variants)
with_ex[f'it_{i}'] = wi
without_ex[f'it_{i}'] = wo

##

# Manual it 1
i += 1
afms = { 
    clone : subset_afm(afm, iteration=i-1, partition=clone) \
    for clone in afm.obs[f'it_{i-1}'].unique() 
}


# One MT-clone, split
a_clone = afms[11] # Take 0
_, a = filter_cells_and_vars(
    a_clone, sample=str(16), filtering='MQuad',
    min_cell_number=0, 
    min_cov_treshold=50, nproc=4, path_=path_tmp
)
a.obs['GBC'].value_counts()
max_n_clones = a.shape[1] if a.shape[1] <= 10 else 10
a.obs[f'it_{i}'] = vireo_wrapper(a, max_n_clones=max_n_clones)
exclusive_variants = { 
    c : rank_clone_variants(a, var=f'it_{i}', group=c, log2_min_perc_ratio=4).index.to_list() \
    for c in a.obs[f'it_{i}'].unique() 
}
wi, wo = test_partitions_variants(exclusive_variants)

pd.crosstab(a.obs[f'it_{1}'], a.obs['GBC']).T

normalized_mutual_info_score(a.obs[f'it_{0}'], a.obs['GBC'])

rank_clone_variants(a, var=f'it_{i}', group=2, min_perc_ratio=1)

##

# Manual it 2
i += 1
afms = { 
    clone : subset_afm(a, iteration=i-1, partition=clone) \
    for clone in a.obs[f'it_{i-1}'].unique() 
}

# One MT-clone, split
a_clone = afms[2] # Take 0
_, a = filter_cells_and_vars(
    a_clone, sample=str(2), filtering='MQuad',
    min_cell_number=0, 
    min_cov_treshold=50, nproc=4, path_=path_tmp
)
_, a = filter_cells_and_vars(
    a_clone, sample=str(2), filtering='MQuad',
    min_cell_number=0, 
    min_cov_treshold=50, nproc=4, path_=path_tmp
)
max_n_clones = a.shape[1] if a.shape[1] <= 20 else 20
a.obs[f'it_{i}'] = vireo_wrapper(a, max_n_clones=10)
exclusive_variants = { 
    c : rank_clone_variants(a, var=f'it_{i}', group=c, min_perc_ratio=1.5).index.to_list() \
    for c in a.obs[f'it_{i}'].unique() 
}
wi, wo = test_partitions_variants(exclusive_variants)

pd.crosstab(a.obs[f'it_{i}'], a.obs['GBC'])

rank_clone_variants(a, var=f'it_{i}', group=2, min_perc_ratio=1)














# Here we go





























# Score
# gt = a.obs['GBC'].astype('str')
# L = []
# ari = custom_ARI(labels, gt)
# nmi = normalized_mutual_info_score(labels, gt)

# pd.crosstab(labels, gt)
