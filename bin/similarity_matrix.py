#!/usr/bin/python

# Calculate cell-cell distance/affinity matrix script

########################################################################

# Parsing CLI args 

# Libraries
import sys 
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='similarity_matrix',
    description=
        '''
        Calculates a cell-to-cell pairwise distance or affinity matrix for all (included) cells in a sample(s).
        '''
)

# Add arguments

# Path_main
my_parser.add_argument(
    '-p', 
    '--path_main', 
    type=str,
    default='..',
    help='The path to the main project directory. Default: .. .'
)

# Sample
my_parser.add_argument(
    '--sample', 
    type=str,
    default='MDA',
    help='Sample to use. Default: MDA. If None, this is done for all cells in all samples found in $path_main/data.'
)

# Filtering
my_parser.add_argument(
    '--filtering', 
    type=str,
    default='miller2022',
    help='Method to filter MT-SNVs. Default: miller2022.'
)

# metric
my_parser.add_argument(
    '--metric', 
    type=str,
    default='euclidean',
    help='Distance metric chosen. Default: euclidean.'
)

# ncores
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=8,
    help='ncores to use for model training. Default: 8.'
)

# min_cell_number
my_parser.add_argument(
    '--min_cell_number', 
    type=int,
    default=0,
    help='Include in the analysis only cells with membership in clones with >= min_cell_number. Default: 0.'
)

# min_cov_treshold
my_parser.add_argument(
    '--min_cov_treshold', 
    type=int,
    default=30,
    help='Include in the analysis only cells MAESTER sites mean coverage > min_cov_treshold. Default: 30.'
)

# Kernel
my_parser.add_argument(
    '--kernel', 
    type=str,
    default=None,
    help='Transform some distance metric with some kernel. Default: None'
)

# Retain nans
my_parser.add_argument(
    '--retain_nans', 
    type=str,
    default='no',
    help='Retain nans. Default no.'
)

# Skip
my_parser.add_argument(
    '--skip', 
    action='store_true',
    help='Skip analysis. Default: False.'
)

# Parse arguments
args = my_parser.parse_args()

path_main = args.path_main
sample = args.sample
filtering = args.filtering
metric = args.metric
ncores = args.ncores 
kernel = args.kernel if args.kernel is not None else 'no_kernel'
retain_nans = args.retain_nans
min_cell_number = args.min_cell_number
min_cov_treshold = args.min_cov_treshold

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    from Cellula._utils import Timer, set_logger
    from MI_TO.preprocessing import *
    from MI_TO.distances import *

    #-----------------------------------------------------------------#

    # Set other paths
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/distances/'
    path_runs = path_main + '/runs/'

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, f'logs_{sample}_{filtering}_{min_cell_number}_{min_cov_treshold}_{metric}_{retain_nans}_{kernel}.txt')

########################################################################

# Main
def main():

    T = Timer()
    T.start()

    # Load data
    t = Timer()
    t.start()

    logger.info(f'Execute pairwise distances calculations: --sample {sample} --filtering {filtering} --metric {metric} --kernel {kernel} --retain_nans {retain_nans} --min_cell_number {min_cell_number} --min_cov_treshold {min_cov_treshold}')

    # Read data
    afm = read_one_sample(path_main, sample=sample)
    ncells0 = afm.shape[0]
    n_all_clones = len(afm.obs['GBC'].unique())

    # Filter 'good quality' cells and variants
    a_cells, a = filter_cells_and_vars(
        afm, 
        filtering=filtering, 
        min_cell_number=min_cell_number, 
        min_cov_treshold=min_cov_treshold, 
        nproc=ncores, 
        path_=path_results
    )
    
    # Format afm for D computation
    if retain_nans == 'no':
        a = nans_as_zeros(a)
        logger.info('Convert nans into zeros...')
    else:
        logger.info('nans mantained in the filtered feature matrix...')
    ncells = a.shape[0]
    n_clones_analyzed = len(a.obs['GBC'].unique())

    logger.info(f'Reading and formatting AFM: total {t.stop()} s.')
    logger.info(f'Total cells and clones in the original QCed sample (only transcriptional and perturb seq QC metrics): {ncells0}; {n_all_clones}.')
    logger.info(f'Total cells, clones and variants in final filtered sample: {ncells}; {n_clones_analyzed}, {a.shape[1]}.')

    # Calculate distance matrix
    t.start()
    logger.info('Begin pairwise distances calculations...')
    D = pair_d(a.X, metric=metric, ncores=ncores, nans=True if retain_nans == 'yes' else False)
    logger.info(f'Finished with distances: {t.stop()} s.')

    # Save as .h5ad adata object
    df_cells = pd.DataFrame({'cell' : a.obs_names}).set_index('cell')
    D = anndata.AnnData(X=D, obs=df_cells, var=df_cells)
    D.uns['distance_matrix'] = {
        'filtering' : filtering, 
        'retain_nans' : retain_nans, 
        'metric' : metric, 
        'kernel' : kernel, 
        'min_cell_number' : min_cell_number, 
        'min_cov_treshold' : min_cov_treshold
    }
    D.write(path_results + f'{sample}_{filtering}_{min_cell_number}_{min_cov_treshold}_{retain_nans}_{metric}_{kernel}.h5ad')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        main()

#######################################################################


