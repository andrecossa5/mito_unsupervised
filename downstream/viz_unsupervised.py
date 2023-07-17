"""
Viz best models and their labels.
"""

# Code
import os
import pickle
from mito_utils.preprocessing import *
from mito_utils.clustering import *
from mito_utils.dimred import *
from mito_utils.embeddings_plots import *
from mito_utils.utils import *
from mito_utils.colors import *
from mito_utils.plotting_base import *
import matplotlib
matplotlib.use('macOSX')


##

# Paths
path_main = '/Users/IEO5505/Desktop/mito_bench/'

# Set paths
path_data = os.path.join(path_main, 'data')
path_output = os.path.join(path_main, 'results', 'unsupervised_clones', 'output')
path_report = os.path.join(path_main, 'results', 'unsupervised_clones', 'reports')
path_viz = '/Users/IEO5505/Desktop/MI_TO/images_DIPA/images'

# Load report
df = pd.read_csv(os.path.join(path_report, 'report_unsupervised.csv'), index_col=0)

# Get top models
df_ = (
    df.loc[:, ['NMI', 'ARI', 'model', 'filtering', 'sample']]
    .assign(
        model_short=lambda x: np.where(x['model'].str.startswith('leiden'), 'leiden', 'vireoSNP'),
        method=lambda x: x['filtering'] + '_' + x['model_short']
    )
    .drop(columns=['model', 'filtering', 'model_short'])
    .groupby('sample')
    .apply(lambda x: x.sort_values('NMI', ascending=False).head(1))
    .reset_index(drop=True)
    .set_index('sample')
)


##


# Master figure!!
fig, axs = plt.subplots(2,4, figsize=(15, 7.8))

# Here we go...
for i, sample in enumerate(['AML_clones', 'MDA_clones', 'MDA_lung', 'MDA_PT']):

    # Get best model: variants, NMI and ARI, GT, labels
    sample_d = df_.loc[sample, :].to_dict()

    # Get unsupervised performance
    with open(os.path.join(path_output, f'out_vireo_{sample}_MQuad_10.pickle'), 'rb') as f:
        d_uns = pickle.load(f)

    # Get variants
    path_ = os.path.join(
        path_data, '../results', 'supervised_clones', 'output',
        f'out_{sample}_MQuad_no_dimred_xgboost_random_10.pickle'
    )
    with open(path_, 'rb') as f:
        d_sup = pickle.load(f)

    d_ = d_sup['trained_models']
    clone = list(d_.keys())[0]
    variants = d_[clone]['variants']

    # Get AFM
    afm = read_one_sample(path_data, sample=sample, with_GBC=True)
    _, a = filter_cells_and_vars(
        afm,
        sample=sample,
        min_cell_number=10, 
        min_cov_treshold=50,
        variants=variants
    )
    a = nans_as_zeros(a)

    # UMAP 
    assert len(d_uns['labels']) == a.shape[0]
    X, _ = reduce_dimensions(a, 'UMAP', n_comps=2, seed=2)
    embs = pd.DataFrame(X, columns=['UMAP1', 'UMAP2'], index=a.obs_names).join(a.obs)
    embs['inference'] = pd.Series(d_uns['labels'], index=a.obs.index)
    embs['GBC'] = embs['GBC'].astype('str')

    # Colors
    colors_gbc, colors_inference = harmonize_colors(embs)
    
    len(colors_inference)
    len(colors_gbc.keys())

    # Axes
    draw_embeddings(embs, cat='GBC', ax=axs[0,i], 
        title=f'{sample}: lentiviral clones',
        axes_kwargs={'legend':False},
        legend_kwargs={'colors':colors_gbc}
    )
    #axs[0,i].axis('off')
    axs[0,i].text(.88, .05, f'n: {int(embs["GBC"].unique().size)}', transform=axs[0,i].transAxes)
    annot_dict = d_uns['df_performance'].T[0].to_dict()
    draw_embeddings(embs, cat='inference', ax=axs[1,i], 
        title=f'{sample}: inference (NMI {annot_dict["NMI"]:.2f}, ARI {annot_dict["ARI"]:.2f})',
        axes_kwargs={'legend':False},
        legend_kwargs={'colors':colors_inference}        
    )
   # axs[1,i].axis('off')
    axs[1,i].text(.88, .05, f'n: {int(embs["inference"].unique().size)-1}', transform=axs[1,i].transAxes)

##

# Viz
fig.tight_layout()
fig.savefig(os.path.join(path_viz, 'unsupervised_UMAPs.png'), dpi=500)
      
##


