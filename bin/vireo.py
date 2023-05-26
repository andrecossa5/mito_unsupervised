#!/usr/bin/python

# Clonal inference with VireoSNP

########################################################################

# Code
import argparse
import os

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
    '--path_data', 
    type=str,
    default='..',
    help='Path input file. Default: .. .'
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

# Chosen 
my_parser.add_argument(
    '--chosen', 
    type=int,
    default=None,
    help='Chosen n of clones, overrides automatic identification via findknee.'
)

# Range n
my_parser.add_argument(
    '--filtering', 
    type=str,
    default='MQuad',
    help='Method to filter variants from the AFM. Default: 2:15.'
)

# min_cov_treshold
my_parser.add_argument(
    '--min_cov_treshold', 
    type=int,
    default=50,
    help='min_cov_treshold.'
)

# min_cell_number
my_parser.add_argument(
    '--min_cell_number', 
    type=int,
    default=0,
    help='min_cell_number treshold.'
)

# min_cell_number
my_parser.add_argument(
    '--p_treshold', 
    type=float,
    default=0.8,
    help='Treshold use to convert clonal assignment to crisp labels.'
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


##

path_data = args.path_data
sample = args.sample
start, stop = args.range.split(':')
range_clones = range(int(start), int(stop))
chosen = args.chosen
filtering = args.filtering
min_cov_treshold = args.min_cov_treshold
min_cell_number = args.min_cell_number
p_treshold = args.p_treshold
ncores = args.ncores

# path_data = '/Users/IEO5505/Desktop/mito_bench/data'
# path_ = os.getcwd()
# sample = 'MDA_PT'
# 
# start, stop = '2:100'.split(':') # 50 stop per PT
# range_clones = range(int(start), int(stop))
# filtering = 'MQuad'
# min_cov_treshold = 50
# min_cell_number = 10
# p_treshold = 0.85
# ncores = 8


##


####################################################################

# Preparing run: import code, set paths

# Code
from kneed import KneeLocator
from sklearn.metrics import normalized_mutual_info_score
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
from mito_utils.embeddings_plots import *
from mito_utils.dimred import *
from mito_utils.clustering import *
from vireoSNP import BinomMixtureVB
import warnings
warnings.filterwarnings("ignore")

####################################################################

# Main
def main():

    # Load data
    T = Timer()
    T.start()
    
    t = Timer()

    # Logging
    path_ = os.getcwd()
    logger = set_logger(path_, f'log_vireo_{sample}_{filtering}_{min_cell_number}.txt')
    
    # Set logger
    logger.info(
        f""" 
        Execute vireoSNP: \n
        --sample {sample} 
        --filtering {filtering} 
        --min_cell_number {min_cell_number} 
        --min_cov_treshold {min_cov_treshold}
        --p_treshold {p_treshold}
        """
    )
    
    ##
    
    # Read data
    t.start()
    afm = read_one_sample(path_data, sample=sample, with_GBC=True)

    # Get variants for MDA_PT
    # path_output = os.path.join(path_data, '../results', 'supervised_clones', 'output')
    # with open(os.path.join(path_output, 'out_MDA_PT_MQuad_no_dimred_xgboost_random_10.pickle'), 'rb') as f:
    #     d = pickle.load(f)
    # variants = d['trained_models']['GGTGACCGCCATCCAATG_vs_rest'] ['variants']

    # Filter cells and vars, create a mutational embedding
    _, a = filter_cells_and_vars(
        afm,
        sample=sample,
        filtering=filtering, 
        min_cell_number=min_cell_number,
        min_cov_treshold=min_cov_treshold,
        #variants=variants,
        nproc=ncores, 
        path_=path_
    )

    # Extract filtered feature matrix, format and reduce with UMAP
    a = nans_as_zeros(a) # For sklearn APIs compatibility
    
    # Get parallel matrices
    AD, DP, _ = get_AD_DP(a, to='csc')
    logger.info(f'SNVs filtering, UMAP, AD/DP matrix preparation {t.stop()}')

    ##

    # Choose k
    t.start()
    logger.info(f'Start inference...')

    _ELBO_mat = []
    for k in range_clones:
        print(f'Clone n: {k}')
        _model = BinomMixtureVB(n_var=AD.shape[0], n_cell=AD.shape[1], n_donor=k)
        _model.fit(AD, DP, min_iter=30, max_iter=500, max_iter_pre=250, n_init=50, random_seed=1234)
        _ELBO_mat.append(_model.ELBO_inits)
    logger.info(f'Finished inference: {t.stop()}')
    
    # Find n_clones
    t.start()
    x = range_clones
    y = np.median(_ELBO_mat, axis=1)
    knee = KneeLocator(x, y).find_knee()[0]
    n_clones = knee
    
    # Refit with optimal n_clones
    _model = BinomMixtureVB(n_var=AD.shape[0], n_cell=AD.shape[1], n_donor=n_clones)
    _model.fit(AD, DP, min_iter=30, n_init=50, max_iter=500, max_iter_pre=250, random_seed=1234)
    
    ##
    
    # Clonal assignment probabilites --> to crisp labels
    clonal_assignment = _model.ID_prob
    df_ass = pd.DataFrame(
        clonal_assignment, 
        index=a.obs_names, 
        columns=range(clonal_assignment.shape[1])
    )

    # Define labels
    labels = []
    for i in range(df_ass.shape[0]):
        cell_ass = df_ass.iloc[i, :]
        try:
            labels.append(np.where(cell_ass>p_treshold)[0][0])
        except:
            labels.append('unassigned')
            
    logger.info(f'Cells assigned to {np.unique(labels).size} clones: {t.stop()}')
        
    # Score
    gt = a.obs['GBC'].astype('str')
    L = []
    ari = custom_ARI(labels, gt)
    nmi = normalized_mutual_info_score(labels, gt)

    L.append({
        'model' : 'vireoSNP',
        'ARI' : ari,
        'NMI' : nmi,
        'sample' : sample,
        'filtering': filtering,
        'dimred' : 'no_dimred',
        'min_cell_number' : 10,
        'min_cov_treshold' : 50,
        'ncells_sample' : gt.size,
        'n_clones_analyzed' : gt.unique().size,
        'n_clones_inferred' : np.unique(labels).size,
        'n_variants' : a.shape[1]
    })
        
    # Write
    df = pd.DataFrame(L)
    d = {'df_performance' : df, 'labels' : pd.Series(labels, index=a.obs_names)}
    with open(
        os.path.join(
            #'/Users/IEO5505/Desktop/mito_bench/results/unsupervised_clones/output',
            f'out_vireo_{sample}_{filtering}_{min_cell_number}.pickle'
            ), 'wb'
        ) as f:
        pickle.dump(d, f)
                   
    # Exit
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

####################################################################

# Run program
if __name__ == "__main__":
    main()
    



