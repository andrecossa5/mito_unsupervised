"""
Trial 2. Take-out-the best simple approach.
"""

# Code
import os
import warnings
from sklearn.metrics import normalized_mutual_info_score
warnings.filterwarnings('ignore')

from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.clustering import *
from mito_utils.iterative import *
from mito_utils.plotting_base import *
from mito_utils.heatmaps_plots import plot_heatmap
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/mito_bench/'
path_data = os.path.join(path_main, 'data') # Here subsampled
path_output = os.path.join(path_main, 'results', 'unsupervised_clones', 'output')
path_viz = os.path.join(path_main, 'results', 'unsupervised_clones', 'visualization')
path_tmp = os.path.join(path_main, 'results', 'unsupervised_clones', 'downstream_files')

##

# AFM first step, no filtering of low numebr of cells clones + variants filtering
afm = read_one_sample(path_data, 'MDA_PT', with_GBC=True)

##

i = 0 
one_dist =  True
final_labels = pd.Series(index=afm.obs_names)
n_variants_final = {}

while one_dist:

    # Filter and cluster
    print(f'Beiginning with {i} iteration...')

    try:

        a_cells, a = filter_cells_and_vars(
            afm if i == 0 else a_cells, # First all cells and vars, then all the good ones.
            filtering='MQuad', 
            min_cov_treshold=50, 
            path_=os.getcwd()
        )

        labels = vireo_wrapper(a)
        a.obs.loc[:, 'g'] = labels

        print(f'n MT-clusters: {labels.unique().size}')

        # Rank vars and find top clone, if any
        n_vars = pd.Series({
            g : rank_clone_variants(
                a, var='g', group=g, rank_by='custom_perc_tresholds',
                min_clone_perc=.9, max_perc_rest=.1
            ).shape[0] \
            for g in a.obs['g'].unique()
        })
        one_dist = np.sum(n_vars>0)>0

        if one_dist:
            top_clone = n_vars.sort_values(ascending=False).index[0]
            n_variants_final[top_clone] = n_vars.loc[top_clone]
            not_top_cells = labels.loc[lambda x: x != top_clone].index
            a_cells, _ = filter_cells_and_vars(
                afm, cells=not_top_cells, variants=a_cells.var_names
            )
            i += 1
            final_labels.loc[labels.loc[lambda x: x == top_clone].index] = str(i)
        else:
            print(f'Finished with {i} iterations. \n')

    except:
        raise ValueError('MQuad or vireo_vrapper ecountered errors. Finishing iterations.')


##


# Process and viz
df = afm.obs.join(
    final_labels.loc[lambda x: x.notna()].to_frame('MT_clones'),
    how='right'
)
df.to_csv(os.path.join(path_tmp, 'trial_2_PT_full.csv'))


##


# Fig   
fig, axs = plt.subplots(2,1,figsize=(8,8))

ari_all = custom_ARI(df['GBC'].astype('str'), df['MT_clones'].astype('str'))
nmi_all = normalized_mutual_info_score(df['GBC'].astype('str'), df['MT_clones'].astype('str'))
df_ = pd.crosstab(df['MT_clones'].astype('str'), df['GBC'].astype('str'))
plot_heatmap(df_, ax=axs[0], y_names_size=5, rank_diagonal=True, label='n cells')
format_ax(title=f'NMI: {nmi_all:.2f}, ARI: {ari_all:.2f}', ax=axs[0], rotx=90)

# > 10 cells
top_clones = df.groupby('GBC').size().loc[lambda x: x>=10].index
df = df.query('GBC in @top_clones')
ari_top = custom_ARI(df['GBC'].astype('str'), df['MT_clones'].astype('str'))
nmi_top = normalized_mutual_info_score(df['GBC'].astype('str'), df['MT_clones'].astype('str'))

df_ = pd.crosstab(df['MT_clones'].astype('str'), df['GBC'].astype('str'))
plot_heatmap(df_, ax=axs[1], y_names_size=5, rank_diagonal=True, label='n cells')
format_ax(title=f'NMI: {nmi_top:.2f}, ARI: {ari_top:.2f}', ax=axs[1], rotx=90)
fig.tight_layout()

fig.savefig(os.path.join(path_output, 'trial_2_PT_full.png'))




