"""
Diagnostics unsupervised MDA_PT: MQuad-vireoSNP
"""

# Code
import os
import re
import pickle
from mito_utils.preprocessing import *
from mito_utils.clustering import *
from mito_utils.dimred import *
from mito_utils.embeddings_plots import *
from mito_utils.utils import *
from mito_utils.colors import *
from mito_utils.plotting_base import *
from mito_utils.heatmaps_plots import *
import matplotlib
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/mito_bench/'
path_data = os.path.join(path_main, 'data')
path_output = os.path.join(path_main, 'results', 'unsupervised_clones', 'output')
path_report = os.path.join(path_main, 'results', 'unsupervised_clones', 'reports')
path_viz = os.path.join(path_main, 'results', 'unsupervised_clones', 'visualization')
path_tmp = os.path.join(path_main, 'results', 'unsupervised_clones', 'downstream_files')

##

# Load report
df = pd.read_csv(os.path.join(path_report, 'report_unsupervised.csv'), index_col=0)

# df.query('sample == "MDA_PT"').sort_values('NMI', ascending=False)

# Pick up labels
with open(os.path.join(path_output, 'out_vireo_MDA_PT_MQuad_10.pickle'), 'rb') as f:
    d = pickle.load(f)

# AFM
afm = read_one_sample(path_data, 'MDA_PT', with_GBC=True)
# Vars
with open(
    os.path.join(path_tmp,
    '/Users/IEO5505/Desktop/mito_bench/results/supervised_clones/downstream_files/MDA_PT_variants.pickle')
    , 'rb') as f:
    variants = pickle.load(f)        # First round variants
# Subset 
a_cells, a = filter_cells_and_vars(
    afm, sample='MDA_PT', variants=variants['MQuad'],
    # filtering='MQuad', path_=path_tmp,
    min_cell_number=10, min_cov_treshold=50, nproc=4
)
a = nans_as_zeros(a)

# Contingency matrix
a.obs['inference'] = pd.Series(d['labels'], index=a.obs_names)
df_ = pd.crosstab(a.obs['inference'], a.obs['GBC'], normalize=1)

# Viz
fig, ax = plt.subplots(figsize=(10,6))
plot_heatmap(df_, ax=ax, annot=False, vmin=.2, vmax=.8,
            rank_diagonal=True, title='MDA_PT: inference vs ground truth')
fig.tight_layout()
fig.savefig(
    os.path.join(
        path_viz, 'MDA_PT_contingency_table.png'
    ),
    dpi=500
)

# Choose some MT_clones: lentiviral families
# np.sum(df_>0.3, axis=1).sort_values() # 11 e 16
cells_11 = a.obs.query('inference == 11').index
cells_16 = a.obs.query('inference == 16').index

a_16 = a[cells_16].copy()
a_16.uns

