#!/usr/bin/python

# Leiden clustering script

########################################################################

# Parsing CLI args 

# Libraries
import sys 
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='leiden_clustering',
    description=
    '''
    Takes each sample AFM and partition the resulting kNN graph across some resolutions
    with the leiden algorithm.
    '''
)

# Add arguments

# Path_main
my_parser.add_argument(
    '-p', 
    '--path_data', 
    type=str,
    default='..',
    help='Path to samples data. Default: .. .'
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
    '--filtering', 
    type=str,
    default='ludwig2019',
    help='Method to filter MT-SNVs. Default: ludwig2019.'
)

# Dimred
my_parser.add_argument(
    '--dimred', 
    type=str,
    default='no_dimred',
    help='Method to reduce the dimension of the SNVs space (top 1000 SNVs selected) by pegasus. Default: no_dimred.'
)

# Dimred
my_parser.add_argument(
    '--n_comps', 
    type=int,
    default=30,
    help='n of components in the dimensionality reduction step. Default: 30.'
)

# GS_mode
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=8,
    help='n of cores used for MQuad filtering and silhouette scoring. Default: 8.'
)

# min_cell_number
my_parser.add_argument(
    '--min_cell_number', 
    type=int,
    default=10,
    help='Include in the analysis only cells with membership in clones with >= min_cell_number. Default: 10.'
)

# min_cov_treshold
my_parser.add_argument(
    '--min_cov_treshold', 
    type=int,
    default=50,
    help='Include in the analysis only cells MAESTER sites mean coverage > min_cov_treshold. Default: 50.'
)

# Resolution range
my_parser.add_argument(
    '--range', 
    type=str,
    default='0.2:1',
    help='Resolution range for leiden clustering. Default: 0.2:1.'
)

# n resolution values
my_parser.add_argument(
    '--n', 
    type=int,
    default=10,
    help='n resolution values to partition the kNN graph with. Default: 10.'
)

# min_cov_treshold
my_parser.add_argument(
    '--k', 
    type=int,
    default=30,
    help='k neighbors. Default: 30.'
)

# skip
my_parser.add_argument(
    '--skip', 
    action='store_true',
    help='Skip analysis. Default: False.'
)

# Parse arguments
args = my_parser.parse_args()

path_data = args.path_data
sample = args.sample
dimred = args.dimred
filtering = args.filtering if dimred == 'no_dimred' else 'pegasus'
min_cell_number = args.min_cell_number
min_cov_treshold = args.min_cov_treshold
n_comps = args.n_comps
res_range = res_range = [ float(x) for x in args.range.split(':') ]
ncores = args.ncores
k = args.k
n = args.n

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import pickle
    from sklearn.metrics import normalized_mutual_info_score
    from mito_utils.utils import *
    from mito_utils.preprocessing import *
    from mito_utils.dimred import *
    from mito_utils.clustering import *
    from mito_utils.kNN import *
    from mito_utils.dimred import *
 
    # Set logger
    path_ = os.getcwd()
    logger = set_logger(
        path_, 
        f'log_leiden_{sample}_{filtering}_{dimred}_{min_cell_number}_{k}.txt'
    )

########################################################################

# Main
def main():

    T = Timer()
    T.start()

    # Load data
    t = Timer()
    t.start()

    logger.info(
        f""" 
        Execute leiden clustering: \n
        --sample {sample} 
        --filtering {filtering} 
        --dimred {dimred}
        --min_cell_number {min_cell_number} 
        --min_cov_treshold {min_cov_treshold}
        """
    )
    
    afm = read_one_sample(path_data, sample=sample, with_GBC=True)
    ncells0 = afm.shape[0]
    n_all_clones = len(afm.obs['GBC'].unique())

    ##

    # Filter 'good quality' cells and variants
    if dimred == 'no_dimred':

        _, a = filter_cells_and_vars(
            afm,
            sample=sample,
            filtering=filtering, 
            min_cell_number=min_cell_number, 
            min_cov_treshold=min_cov_treshold, 
            nproc=ncores, 
            path_=path
        )

        # Extract X
        a = nans_as_zeros(a) # For sklearn APIs compatibility
        cells = a.obs_names
        variants = a.var_names
        ncells = cells.size
        n_clones_analyzed = len(a.obs['GBC'].unique())
        X = a.X

    else:

        _, a = filter_cells_and_vars(
            afm,
            sample=sample,
            filtering=filtering, 
            min_cell_number=min_cell_number, 
            min_cov_treshold=min_cov_treshold, 
            nproc=ncores,
            path_=path_
        )

        # Extract X
        a = nans_as_zeros(a) # For sklearn APIs compatibility
        cells = a.obs_names
        variants = a.var_names
        ncells = cells.size
        n_clones_analyzed = len(a.obs['GBC'].unique())
        X, _ = reduce_dimensions(a, method=dimred, n_comps=min(n_comps, ncells), sqrt=False)

    # Here we go
    logger.info(f'Finished cells and variant selection {t.stop()}')
    logger.info(f'Execute leiden clustering...')
    
    # kNN
    t.start()
    g = kNN_graph(X, k=k)
    conn = g[1]
    logger.info(f'kNN construction: {t.stop()}')

    # Partition and scoring
    gt = a.obs['GBC'].astype('str')
    labels_d = {}
    L = []
    for i, r in enumerate(np.linspace(res_range[0], res_range[1], n)):

        t.start()
        labels = leiden_clustering(conn, res=r)
        ari = custom_ARI(labels, gt)
        nmi = normalized_mutual_info_score(labels, gt)
        
        L.append({
            'model' : f'leiden_k_{k}_res_{r}',
            'ARI' : ari,
            'NMI' : nmi,
            'sample' : sample,
            'filtering': filtering,
            'dimred' : dimred,
            'min_cell_number' : min_cell_number,
            'min_cov_treshold' : min_cov_treshold,
            'ncells_sample' : gt.size,
            'n_clones_analyzed' : gt.unique().size,
            'n_clones_inferred' : np.unique(labels).size,
            'n_variants' : a.shape[1] if dimred != 'no_dimred' else 'reduced_dim',
        })
        labels_d[r] = labels
        logger.info(f'Sol {k}_{r}: {t.stop()}')

    # Write 
    df = pd.DataFrame(L).sort_values('NMI')
    d = {
        'df_performance' : df,
        'labels' : labels_d
    }
    path_results = f'out_leiden_{sample}_{filtering}_{dimred}_{min_cell_number}_{k}.pickle'
    with open(path_results, 'wb') as f:
        pickle.dump(d, f)
    
    #-----------------------------------------------------------------#

    # Exec
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')


########################################################################   

if __name__ == "__main__":
    main()