"""
Script to compute leiden clustering visualizations.
"""

import pickle
from Cellula._utils import Timer, set_logger
from Cellula.preprocessing._metrics import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._colors import *
from MI_TO.utils import *
from MI_TO.preprocessing import *
from MI_TO.kNN import *
from MI_TO.spectral_clustering import *
from MI_TO.dimensionality_reduction import *
matplotlib.use('macOSX')

# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/'
path_data = path_main + '/data/'
path_distances = path_main + '/results_and_plots/distances/'
path_clones = path_main + '/results_and_plots/clones_classification/'
path_results = path_main + '/results_and_plots/leiden_clustering/'
path_class_performance = path_main + '/results_and_plots/classification_performance/'
path_runs = path_main + '/runs/'

# Read dictionaries
with open(path_clones + 'top3.pkl', 'rb') as f:
    top3 = pickle.load(f)

with open(path_results + 'solutions.pkl', 'rb') as f:
    solutions = pickle.load(f)


with open(path_results + 'connectivities.pkl', 'rb') as f:
    connectivities = pickle.load(f)

# Extract all info
analysis_l = []
metric_l = []
ARI_l = []
sample_l = []
k_l = []
res_l = []
id_l = []
with_nans_l = []

for sample in solutions:
    d_ = solutions[sample]
    for x in d_:
        id_l.append(x)
        similarity, k, res = x.split('|')
        analysis_l.append('_'.join(similarity.split('_')[:4]))
        sample_l.append(similarity.split('_')[0])
        metric_l.append(similarity.split('_')[5])
        with_nans_l.append(similarity.split('_')[4])
        labels, cells, ARI = d_[x]
        ARI_l.append(ARI)
        k_l.append(k)
        res_l.append(str(round(float(res), 2)))

df = pd.DataFrame(
    {   
        'id' : id_l,
        'analysis' : analysis_l,
        'sample' : sample_l,
        'metric' : metric_l,
        'with_nans' : with_nans_l,
        'k' : k_l,
        'res' : res_l,
        'ARI' : ARI_l,
    }
).set_index('id')


##


###############Summarise leiden clustering concordance (ARI)

fig, axs = plt.subplots(2,2,figsize=(10,8), sharey=True)

params = {   
    'showcaps' : True,
    'fliersize': 0,
    'boxprops' : {'edgecolor': 'black', 'linewidth': 0.3}, 
    'medianprops': {"color": "black", "linewidth": 1},
    'whiskerprops':{"color": "black", "linewidth": 1}
}

box(df, 'metric', 'ARI', c='#E9E7E7', params=params, ax=axs[0,0])
strip(df, 'metric', 'ARI', c='#000066', s=4, ax=axs[0,0])
format_ax(df, ylabel='ARI', title='ARI by metric', ax=axs[0,0])

box(df, 'k', 'ARI', c='#E9E7E7', params=params, ax=axs[0,1])
strip(df, 'k', 'ARI', c='#000066', s=4, ax=axs[0,1])
format_ax(df, title='ARI by n NN neighbors', ax=axs[0,1])

box(df, 'res', 'ARI', c='#E9E7E7', params=params, ax=axs[1,0,])
strip(df, 'res', 'ARI', c='#000066', s=4, ax=axs[1,0])
format_ax(df, ylabel='ARI', title='ARI by leiden clustering resolution', ax=axs[1,0])

box(df, 'sample', 'ARI', c='#E9E7E7', params=params, ax=axs[1,1])
strip(df, 'sample', 'ARI', c='#000066', s=4, ax=axs[1,1])
format_ax(df, title='ARI by sample',ax=axs[1,1])

fig.tight_layout()
fig.savefig(path_results + 'summary_ARI.pdf')

###############


##


# Compute top concordant solutions UMAPs embeddings
top_runs_per_sample = df.sort_values(by='ARI', ascending=False).groupby('sample').head(1)

############### MDA umap

# Prep umap
X, conn, cells, true_clones, labels, ARI, d_run = prep_things_for_umap(top_runs_per_sample, 1, solutions, connectivities, path_main=path_main)

# Calculate UMAP
X_umap = umap_from_X_conn(X, conn, metric=d_run['metric'])
df_ = pd.DataFrame(X_umap, columns=['UMAP1', 'UMAP2'], index=cells)
df_['true_clones'] = true_clones
df_['leiden'] = labels
df_['leiden'] = df_['leiden'].astype(str)

df_

# Colors
with open(path_main + 'data/MDA_clones_colors.pkl', 'rb') as f:
    true_colors = pickle.load(f)
true_colors = { k:v for k,v in true_colors.items() if k in df_['true_clones'].unique() }
leiden_colors = create_palette(df_, 'leiden', sc.pl.palettes.vega_10)

# Viz
fig, axs = plt.subplots(1,2, figsize=(10,5), constrained_layout=True)

scatter(df_, 'UMAP1', 'UMAP2', by='true_clones', c=true_colors, ax=axs[0], s=7)
format_ax(df_, ax=axs[0], title='Ground truth')
handles = create_handles(true_colors.keys(), colors=true_colors.values())
axs[0].legend(handles, true_colors.keys(), loc='center', 
    bbox_to_anchor=(0.25, 0.18), ncol=1, frameon=False, title='Perturb-seq clones', fontsize=9
)
axs[0].axis('off')

scatter(df_, 'UMAP1', 'UMAP2', by='leiden', c=leiden_colors, ax=axs[1], s=7)
format_ax(df_, ax=axs[1], title='Leiden clusters')
handles = create_handles(leiden_colors.keys(), colors=leiden_colors.values())
axs[1].legend(handles, leiden_colors.keys(), loc='center', 
    bbox_to_anchor=(0.25, 0.18), ncol=1, frameon=False, title='Leiden clusters', fontsize=9
)
axs[1].axis('off')

fig.suptitle(f'MDA sample: top solution ARI {ARI:.2f}')
fig.savefig(path_results + 'MDA_top.pdf')

###############


##


############### AML umap

# Prep umap
X, conn, cells, true_clones, labels, ARI, d_run = prep_things_for_umap(top_runs_per_sample, 0, solutions, connectivities, path_main=path_main)

# Calculate UMAP
X_umap = umap_from_X_conn(X, conn, metric=d_run['metric'])
df_ = pd.DataFrame(X_umap, columns=['UMAP1', 'UMAP2'], index=cells)
df_['true_clones'] = true_clones
df_['leiden'] = labels
df_['leiden'] = df_['leiden'].astype(str)

df_

# Colors
with open(path_main + 'data/AML_clones_colors.pkl', 'rb') as f:
    true_colors = pickle.load(f)
true_colors = { k:v for k,v in true_colors.items() if k in df_['true_clones'].unique() }
leiden_colors = create_palette(df_, 'leiden', sc.pl.palettes.vega_10)

# Viz
fig, axs = plt.subplots(1,2, figsize=(10,5), constrained_layout=True)

scatter(df_, 'UMAP1', 'UMAP2', by='true_clones', c=true_colors, ax=axs[0], s=7)
format_ax(df_, ax=axs[0], title='Ground truth')
handles = create_handles(true_colors.keys(), colors=true_colors.values())
axs[0].legend(handles, true_colors.keys(), loc='center', 
    bbox_to_anchor=(0.25, 0.18), ncol=1, frameon=False, title='Perturb-seq clones', fontsize=9
)
axs[0].axis('off')

scatter(df_, 'UMAP1', 'UMAP2', by='leiden', c=leiden_colors, ax=axs[1], s=7)
format_ax(df_, ax=axs[1], title='Leiden clusters')
handles = create_handles(leiden_colors.keys(), colors=leiden_colors.values())
axs[1].legend(handles, leiden_colors.keys(), loc='center', 
    bbox_to_anchor=(0.25, 0.18), ncol=1, frameon=False, title='Leiden clusters', fontsize=9
)
axs[1].axis('off')

fig.suptitle(f'AML sample: top solution ARI {ARI:.2f}')
fig.savefig(path_results + 'AML_top.pdf')

###############


##


############### PDX umap

# Prep umap
X, conn, cells, true_clones, labels, ARI, d_run = prep_things_for_umap(top_runs_per_sample, 2, solutions, connectivities, path_main=path_main)

# Calculate UMAP
X_umap = umap_from_X_conn(X, conn, metric=d_run['metric'])
df_ = pd.DataFrame(X_umap, columns=['UMAP1', 'UMAP2'], index=cells)
df_['true_clones'] = true_clones
df_['leiden'] = labels
df_['leiden'] = df_['leiden'].astype(str)

df_

# Colors
with open(path_main + 'data/PDX_clones_colors.pkl', 'rb') as f:
    true_colors = pickle.load(f)
true_colors = { k:v for k,v in true_colors.items() if k in df_['true_clones'].unique() }
leiden_colors = create_palette(df_, 'leiden', sc.pl.palettes.vega_10)

# Viz
fig, axs = plt.subplots(1,2, figsize=(10,5), constrained_layout=True)

scatter(df_, 'UMAP1', 'UMAP2', by='true_clones', c=true_colors, ax=axs[0], s=7)
format_ax(df_, ax=axs[0], title='Ground truth')
handles = create_handles(true_colors.keys(), colors=true_colors.values(), size=6)
axs[0].legend(handles, true_colors.keys(), loc='center', 
    bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False, title='Perturb-seq clones', fontsize=6
)
axs[0].axis('off')

scatter(df_, 'UMAP1', 'UMAP2', by='leiden', c=leiden_colors, ax=axs[1], s=7)
format_ax(df_, ax=axs[1], title='Leiden clusters')
handles = create_handles(leiden_colors.keys(), colors=leiden_colors.values(), size=6)
axs[1].legend(handles, leiden_colors.keys(), loc='center', 
    bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False, title='Leiden clusters', fontsize=6
)
axs[1].axis('off')

fig.suptitle(f'PDX sample: top solution ARI {ARI:.2f}')
fig.savefig(path_results + 'PDX_top.pdf')

################