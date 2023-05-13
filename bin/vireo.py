#!/usr/bin/python

# Clonal inference with VireoSNP

########################################################################

# Code
import argparse
import os
from sklearn.metrics import normalized_mutual_info_score
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.dimred import *
from mito_utils.clustering import *
from vireoSNP import BinomMixtureVB
from vireoSNP.plot import heat_matrix
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


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

# Filter
my_parser.add_argument(
    '--sample', 
    type=str,
    default='MDA_clones',
    help='Sample to use. Default: MDA_clones.'
)

# Filter
my_parser.add_argument(
    '--chosen', 
    type=int,
    default=None,
    help='Chosen n of clones.'
)

# Parse arguments
args = my_parser.parse_args()

path_main = args.path_main
sample = args.sample
chosen = args.chosen

####################################################################

# Paths 
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'vireoSNP')

# Default args here
filtering = 'MQuad'
min_cov_treshold = 50
ncores = 10
p_treshold = 0.8

##

####################################################################

# Main
def main():

    # Move to vireoSNP directory
    T = Timer()
    T.start()

    make_folder(path_results, sample, overwrite=True)
    logger = set_logger(os.path.join(path_results, sample), 'log.txt')
    
    # Read data
    afm = read_one_sample(path_data, sample=sample, with_GBC=False)
    # afm = read_one_sample(path_data, sample=sample)
    # ncells0 = afm.shape[0]
    # n_all_clones = len(afm.obs['GBC'].unique())

    # Filter cells and vars, create a mutational embedding
    _, a = filter_cells_and_vars(
        afm,
        # blacklist=blacklist,
        sample=sample,
        filtering=filtering, 
        min_cov_treshold=min_cov_treshold, 
        nproc=ncores, 
        path_=os.path.join(path_results, sample)
    )

    # Extract filtered feature matrix, format and reduce with UMAP
    a = nans_as_zeros(a) # For sklearn APIs compatibility
    ncells = a.shape[0]

    # UMAP MT-SNVs
    embs, _ = reduce_dimensions(a, method='UMAP', n_comps=min(30, a.shape[1]), sqrt=False)
    (
        pd.DataFrame(embs, index=a.obs_names, columns=['UMAP1', 'UMAP2'])
        .to_csv(os.path.join(path_results, sample, 'MT_SNVs_umap.csv'))
    )

    # Get parallel matrices
    AD, DP, _ = get_AD_DP(a, to='csc')

    ##

    # Choose k
    if chosen is None:

        n_init = 50
        n_clone_list = np.arange(2, 15)

        _ELBO_mat = []
        for k in n_clone_list:
            _model = BinomMixtureVB(n_var=AD.shape[0], n_cell=AD.shape[1], n_donor=k)
            _model.fit(AD, DP, min_iter=30, n_init=n_init)
            _ELBO_mat.append(_model.ELBO_inits)

        # Show ELBO trend
        fig = plt.figure(figsize=(10, 5))
        plt.plot(np.arange(1, len(n_clone_list)+1), np.max(_ELBO_mat, axis=1))
        plt.boxplot(_ELBO_mat)
        plt.xticks(np.arange(1, len(n_clone_list)+1), n_clone_list)
        plt.ylabel("ELBO")
        plt.xlabel("n_clones")

        # Save
        fig.savefig(os.path.join(path_results, sample, 'ELBO.png'))

    ##

    # Identify clones
    if chosen is not None:

        # Last fit
        _model = BinomMixtureVB(n_var=AD.shape[0], n_cell=AD.shape[1], n_donor=chosen)
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
        # ground_truth = a.obs['GBC'].cat.codes.values
        labels = []
        for i in range(df_ass.shape[0]):
            cell_ass = df_ass.iloc[i, :]
            try:
                labels.append(np.where(cell_ass>p_treshold)[0][0])
            except:
                labels.append('unassigned')
        # 
        # 
        # custom_ARI(ground_truth, labels)
        # normalized_mutual_info_score(ground_truth, labels)

        # np.unique(ground_truth, return_counts=True)

        # Save labels
        print(np.unique(labels, return_counts=True))
        (
            pd.Series(labels, index=a.obs_names, columns=['vireoSNP_labels'])
            .to_csv(os.path.join(path_results, sample, 'lables.csv'))
        )

        ##

        # # Model fitting diagnostics 1
        # import matplotlib
        # matplotlib.use('macOSX')
        # 
        # fig = plt.figure(figsize=(11, 4))
        # plt.subplot(1, 2, 1)
        # plt.hist(_model.ELBO_inits)
        # plt.ylabel("Frequency")
        # plt.xlabel("ELBO in multiple initializations")
        # 
        # plt.subplot(1, 2, 2)
        # plt.plot(_model.ELBO_iters)
        # plt.xlabel("Iterations")
        # plt.ylabel("ELBO in a single initialization")
        # 
        # plt.tight_layout()
        # plt.show()

        # Visualize output
        raw_col = cm.get_cmap('pink_r', 200)
        new_col = np.vstack((raw_col(np.linspace(0, 0.7, 10)),
                             raw_col(np.linspace(0.7, 1, 90))))
        segpink = ListedColormap(new_col, name='segpink')

        # Clonal assignment visualization
        fig = plt.figure(figsize=(7, 4))
        plt.subplot(1, 2, 1)
        im = heat_matrix(_model.ID_prob, cmap="Blues", alpha=0.8,
                         display_value=False, row_sort=True)
        plt.colorbar(im, fraction=0.046, pad=0.04)
        plt.title("Assignment probability")
        plt.xlabel("Clone")
        plt.ylabel("%d cells" %(_model.n_cell))
        plt.xticks(range(_model.n_donor))

        plt.subplot(1, 2, 2)
        im = heat_matrix(_model.beta_mu, cmap=segpink, alpha=0.8,
                         display_value=False, row_sort=True)
        plt.colorbar(im, fraction=0.046, pad=0.04)
        plt.title("Mean allelic ratio")
        plt.xlabel("Clone")
        plt.ylabel("%d SNPs" %(_model.n_var))
        plt.xticks(range(_model.n_donor))

        # Save
        plt.tight_layout()
        fig.savefig(os.path.join(path_results, sample, 'heatmaps.png'))

        ##
        
        logger.info(f'Execution was completed successfully in total {T.stop()} s.')

####################################################################

# Run program
if __name__ == "__main__":
    main()