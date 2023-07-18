"""
Trial 2. Take-out-the best simple approach.
"""

# Code
import os
import sys
import warnings
from joblib import cpu_count
from sklearn.metrics import normalized_mutual_info_score
warnings.filterwarnings('ignore')

from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.clustering import *
from mito_utils.iterative import *
from mito_utils.plotting_base import *
from mito_utils.heatmaps_plots import plot_heatmap
# matplotlib.use('macOSX')


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

# Set n for subfiltering
n = int(sys.argv[1]) if sys.argv[1] != 'None' else None


##


def main():

    i = 0 
    one_dist =  True
    final_labels = pd.Series(index=afm.obs_names)
    n_variants_final = {}

    while one_dist:

        # Filter and cluster
        print(f'Beginning {i} iteration...')

        try:

            a_cells, a = filter_cells_and_vars(
                afm if i == 0 else a_cells, # First all cells and vars, then all the good ones.
                filtering='MQuad', 
                min_cov_treshold=50,
                nproc=cpu_count(),
                path_=os.getcwd(),
                n=n
            )
            if i == 0:
                good_quality_cells = a_cells.obs_names

            labels = vireo_wrapper(a)
            a.obs.loc[:, 'g'] = labels

            print(f'n MT-clusters: {labels.unique().size}')

            # Rank vars and find top clone, if any
            t = .7
            while t>.5:
                n_vars = pd.Series({
                    g : rank_clone_variants(
                        a, var='g', group=g, rank_by='custom_perc_tresholds',
                        min_clone_perc=t, max_perc_rest=.1
                    ).shape[0] \
                    for g in a.obs['g'].unique()
                })
                one_dist = np.sum(n_vars>0)>0
                if one_dist:
                    break
                else:
                    t -= .05

            if one_dist:
                top_clone = n_vars.sort_values(ascending=False).index[0]
                n_variants_final[i] = n_vars.loc[top_clone]
                not_top_cells = labels.loc[lambda x: x != top_clone].index
                a_cells, _ = filter_cells_and_vars(
                    afm, cells=not_top_cells, variants=a_cells.var_names
                )
                i += 1
                final_labels.loc[labels.loc[lambda x: x == top_clone].index] = str(i)
            else:
                print(f'Finished with {i} iterations. \n')

        except:
            print('MQuad or vireo_vrapper ecountered errors. Finishing iterations.')
            break


    print(final_labels.value_counts())
    # to_recluster = pd.Series(n_variants_final).loc[lambda x:x>1].index.astype('str')
    afm.obs = afm.obs.join(final_labels.to_frame('g'))

    # g = to_recluster[2]
    # afm.obs['GBC'] = afm.obs['GBC'].astype('str')
    # subset = afm.obs.query('g == @g').index
    # 
    # a_subset = filter_cells_and_vars(afm, cells=subset)[0]
    # a_cells, a = filter_cells_and_vars(
    #     a_subset,
    #     filtering='MQuad', 
    #     min_cov_treshold=50, 
    #     path_=os.getcwd()
    # )
    # labels = vireo_wrapper(a)
    # a.obs.loc[:, 'g'] = labels
    # 
    # print(f'n MT-clusters: {labels.unique().size}')
    # 
    # a.obs.groupby('GBC').size()
    # a.obs.groupby('g').size()
    # rank_clone_variants(
    #     a, var='g', group='1', rank_by='custom_perc_tresholds', filter_vars=False
    #     # min_clone_perc=.2, max_perc_rest=.1
    # )
    # 

    # Rank vars and find top clone, if any
    # t = .7
    # while t>.5:
    #     n_vars = pd.Series({
    #         g : rank_clone_variants(
    #             a, var='g', group=g, rank_by='custom_perc_tresholds',
    #             min_clone_perc=t, max_perc_rest=.1
    #         ).shape[0] \
    #         for g in a.obs['g'].unique()
    #     })
    #     one_dist = np.sum(n_vars>0)>0
    #     if one_dist:
    #         break
    #     else:
    #         t -= .05
    # 
    # if one_dist:
    #     top_clone = n_vars.sort_values(ascending=False).index[0]
    #     n_variants_final[i] = n_vars.loc[top_clone]
    #     not_top_cells = labels.loc[lambda x: x != top_clone].index
    #     a_cells, _ = filter_cells_and_vars(
    #         afm, cells=not_top_cells, variants=a_cells.var_names
    #     )
    #     i += 1
    #     final_labels.loc[labels.loc[lambda x: x == top_clone].index] = str(i)
    # else:
    #     print(f'Finished with {i} iterations. \n')


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
    perc_assigned = df.shape[0] / good_quality_cells.size
    n_clones_gt = afm.obs.loc[good_quality_cells]['GBC'].unique().size
    n_clones_recovered = df['GBC'].unique().size
    t = f'''NMI: {nmi_all:.2f}, ARI: {ari_all:.2f} 
        All clones: {n_clones_recovered}/{n_clones_gt} clones recovered, {perc_assigned*100:.2f}% cells'''
    format_ax(title=t, ax=axs[0], rotx=90)
    print(t)

    # > 10 cells
    top_clones = afm.obs.loc[good_quality_cells].groupby('GBC').size().loc[lambda x: x>=10].index
    n_top_all = afm.obs.loc[good_quality_cells].query('GBC in @top_clones').shape[0]

    df = df.query('GBC in @top_clones')
    ari_top = custom_ARI(df['GBC'].astype('str'), df['MT_clones'].astype('str'))
    nmi_top = normalized_mutual_info_score(df['GBC'].astype('str'), df['MT_clones'].astype('str'))

    df_ = pd.crosstab(df['MT_clones'].astype('str'), df['GBC'].astype('str'))
    plot_heatmap(df_, ax=axs[1], y_names_size=5, rank_diagonal=True, label='n cells')
    n_clones_gt = (
        afm.obs.loc[good_quality_cells]['GBC']
        .value_counts()
        .loc[lambda x: x>=10]  
        .unique().size
    )
    n_clones_recovered = df['GBC'].unique().size
    perc_assigned = (df.shape[0] / n_top_all) * 100
    t = f'''NMI: {nmi_top:.2f}, ARI: {ari_top:.2f} 
            Clones >=10 cells: {n_clones_recovered}/{n_clones_gt} clones recovered, {perc_assigned:.2f}% cells'''
    format_ax(title=t, ax=axs[1], rotx=90)
    print(t)

    fig.tight_layout()
    fig.savefig(os.path.join(path_viz, f'trial_2_PT_subsampled_{n}.png'), dpi=500)


    ##


######################################################

# Run
if __name__ == '__main__':
    main()

