#!/usr/bin/python

# Clonal inference with VireoSNP

########################################################################

# Code
import argparse
import os

from kneed import KneeLocator
from sklearn.metrics import normalized_mutual_info_score
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
from mito_utils.embeddings_plots import *
from mito_utils.dimred import *
from mito_utils.clustering import *
from vireoSNP import BinomMixtureVB


##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='vireoSNP',
    description='Wrapper script within vireoSNP vignette.'
)

# Add arguments

# Path_main
my_parser.add_argument(
    '-p', 
    '--path_main', 
    type=str,
    default='..',
    help='Path to project dir. Default: .. .'
)

# Sample
my_parser.add_argument(
    '--sample', 
    type=str,
    default='MDA_clones',
    help='Sample to use. Default: MDA_clones.'
)

# Range n
my_parser.add_argument(
    '--range', 
    type=str,
    default='2:15',
    help='Range of n_clones to search from. Default: 2:15.'
)

# GBC
my_parser.add_argument(
    '--GBC', 
    action='store_true',
    help='Read and use GBC information. Default: False.'
)

# Chosen 
my_parser.add_argument(
    '--chosen', 
    type=int,
    default=None,
    help='Chosen n of clones, overrides automatic identification via findknee.'
)

# ncores
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=None,
    help='n cores to use.'
)

# Parse arguments
args = my_parser.parse_args()

# path_main = '/Users/IEO5505/Desktop/mito_bench/'
# sample = 'MDA_clones'
# ncores = 4
# chosen = None
# range_clones = np.arange(2, 15)
# GBC = False
# filtering = 'MQuad'
# min_cov_treshold = 50
# dimred
# p_treshold = 0.8

path_main = args.path_main
sample = args.sample
ncores = args.ncores
chosen = args.chosen
start, stop = args.range.split(':')
range_clones = np.arange(int(start), int(stop))
GBC =  args.GBC

# Default
filtering = 'MQuad'
min_cov_treshold = 50
dimred = 'no_dimred'
p_treshold = 0.8

####################################################################

# Paths 
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'unsupervised_clones', 'vireo')

####################################################################

# Main
def main():

    # Move to vireoSNP directory
    T = Timer()
    T.start()

    make_folder(path_results, sample, overwrite=True)
    logger = set_logger(os.path.join(path_results, sample), 'log.txt')
    
    # Read data
    if GBC:
        afm = read_one_sample(path_data, sample=sample, with_GBC=True)
        n_all_clones = len(afm.obs['GBC'].unique())
    else:
        afm = read_one_sample(path_data, sample=sample, with_GBC=False)

    # Filter cells and vars, create a mutational embedding
    ncells0 = afm.shape[0]
    _, a = filter_cells_and_vars(
        afm,
        # blacklist=blacklist,
        sample=sample,
        filtering=filtering, 
        min_cell_number=0 if not GBC else 10,
        min_cov_treshold=min_cov_treshold,
        nproc=ncores, 
        path_=os.path.join(path_results, sample)
    )

    # Extract filtered feature matrix, format and reduce with UMAP
    a = nans_as_zeros(a) # For sklearn APIs compatibility
    ncells = a.shape[0]

    # UMAP MT-SNVs
    embs, _ = reduce_dimensions(a, method='UMAP', n_comps=2, sqrt=False)
    embs = pd.DataFrame(embs, index=a.obs_names, columns=['UMAP1', 'UMAP2'])
    embs.to_csv(os.path.join(path_results, sample, 'MT_SNVs_umap.csv'))

    # Get parallel matrices
    AD, DP, _ = get_AD_DP(a, to='csc')

    ##

    # Choose k
    if chosen is None:

        _ELBO_mat = []
        for k in range_clones:
            _model = BinomMixtureVB(n_var=AD.shape[0], n_cell=AD.shape[1], n_donor=k)
            _model.fit(AD, DP, min_iter=30, max_iter=500, max_iter_pre=250, n_init=50, random_seed=1234)
            _ELBO_mat.append(_model.ELBO_inits)
        
        # Find knee
        x = range_clones
        y = np.median(_ELBO_mat, axis=1)
        knee = KneeLocator(x, y).find_knee()[0]
        n_clones = knee

        # Show ELBO trend with knee
        df_ = (
            pd.DataFrame(np.array(_ELBO_mat), index=range_clones)
            .reset_index()
            .rename(columns={'index':'n_clones'})
            .melt(id_vars='n_clones', var_name='run', value_name='ELBO')
        )    
        df_['n_clones'] = df_['n_clones'].astype(str)
         
        # Fig
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(df_.groupby('n_clones').agg('median')['ELBO'].sort_values(), 'o--')
        ax.vlines(str(knee), df_['ELBO'].min(), df_['ELBO'].max())
        box(df_, 'n_clones', 'ELBO', c='#E9E7E7', ax=ax)
        format_ax(ax, xlabel='n clones', ylabel='ELBO')
        ax.text(.8, .1, f'knee: {knee}', transform=ax.transAxes)
        ax.spines[['right', 'top']].set_visible(False)

        # Save
        fig.savefig(os.path.join(path_results, sample, 'ELBO.png'))

    ##

    # Identify clones
    _model = BinomMixtureVB(n_var=AD.shape[0], n_cell=AD.shape[1], n_donor=n_clones)
    _model.fit(AD, DP, min_iter=30, n_init=50)

    ##
    
    # Clonal assignment probabilites --> to crisp labels
    clonal_assignment = _model.ID_prob
    df_ass = pd.DataFrame(
        clonal_assignment, 
        index=a.obs_names, 
        columns=range(clonal_assignment.shape[1])
    )
    df_ass.to_csv(os.path.join(path_results, sample, 'cell_assignments.csv'))

    # Define labels
    labels = []
    for i in range(df_ass.shape[0]):
        cell_ass = df_ass.iloc[i, :]
        try:
            labels.append(np.where(cell_ass>p_treshold)[0][0])
        except:
            labels.append('unassigned')
    
    # Score if necessary      
    if GBC:
        ground_truth = a.obs['GBC'].astype('str')
        ARI = custom_ARI(ground_truth, labels)
        NMI = normalized_mutual_info_score(ground_truth, labels)

        # Save labels
        print(np.unique(labels, return_counts=True))
        labels = pd.Series(labels, index=a.obs_names)
        labels.to_csv(os.path.join(path_results, sample, 'labels.csv'))
               
    # Exit
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

####################################################################

# Run program
if __name__ == "__main__":
    main()
    



