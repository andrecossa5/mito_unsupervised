#!/usr/bin/python

# Take-one-off scripts

########################################################################

# Parsing CLI args 

# Libraries
import sys 
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='take_one_off_trial',
    description=
    '''
    Take-one-off iterative trial.
    '''
)

# Add arguments

# Path_main
my_parser.add_argument(
    '-p', 
    '--path_main', 
    type=str,
    default='..',
    help='Path to data folder. Default: .. .'
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
    default='MQuad',
    help='Method to filter MT-SNVs. Default: MQuad.'
)

# Dimred
my_parser.add_argument(
    '--n_muts', 
    type=int,
    default=None,
    help='n muts in the prefiltering to MQuad. Default: None, all of them.'
)

# GS_mode
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=8,
    help='n of cores used for MQuad filtering. Default: 8.'
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

# N max muts
my_parser.add_argument(
    '--n_max_mut', 
    action='store_true',
    help='VireoSNP clustering with max k --> n_max_mut. Default: False.'
)

# skip
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
min_cell_number = args.min_cell_number
min_cov_treshold = args.min_cov_treshold
n_muts = args.n_muts
ncores = args.ncores
n_max_mut = args.n_max_mut

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import os
    import warnings
    warnings.filterwarnings('ignore')

    from mito_utils.preprocessing import *
    from mito_utils.utils import *
    from mito_utils.preprocessing import *
    from mito_utils.clustering import *
    from mito_utils.it_diagnostics import *
    from mito_utils.it_iterations import *
    from mito_utils.plotting_base import *
    from mito_utils.it_plots import *
    
    # Set logger
    path_ = os.getcwd()
    make_folder(path_, 'logs', overwrite=False)
    path_ = os.path.join(path_, 'logs')
    logger = set_logger(
        path_, 
        f'logs_take_one_off_{sample}_{filtering}_{n_muts}.txt'
    )

########################################################################

# Main
def main():

    T = Timer()
    T.start()

    t = Timer()

    # Load data
    logger.info(
        f""" 
        Execute iterative clustering: \n
        --sample {sample} 
        --filtering {filtering} 
        --min_cell_number {min_cell_number} 
        --min_cov_treshold {min_cov_treshold}
        """
    )

    ##

    # Set paths
    path_data = os.path.join(path_main, 'data') 
    path_output = os.path.join(path_main, 'results', 'unsupervised_clones', 'output')
    path_viz = os.path.join(path_main, 'results', 'unsupervised_clones', 'visualization')
    path_tmp = os.path.join(path_main, 'results', 'unsupervised_clones', 'downstream_files')

    ##

    # AFM first step, no filtering of low numebr of cells clones + variants filtering
    afm = read_one_sample(path_data, sample, with_GBC=True)

    i = 0 
    one_dist =  True
    final_labels = pd.Series(index=afm.obs_names)
    n_variants_final = {}

    # Here we go
    while one_dist:

        # Filter and cluster
        t.start()
        logger.info(f'Beginning {i} iteration...')

        try:
            a_cells, a = filter_cells_and_vars(
                afm if i == 0 else a_cells, # First all cells and vars, then all the good ones.
                filtering=filtering, 
                min_cov_treshold=50,
                nproc=ncores,
                path_=os.getcwd(),
                n=n_muts
            )
            if i == 0:
                good_quality_cells = a_cells.obs_names

            labels = vireo_wrapper(a, n_max_mut=n_max_mut)
            a.obs.loc[:, 'g'] = labels

            logger.info(f'n MT-clusters: {labels.unique().size}')

            # Rank vars and find top clone, if any
            n_vars = find_exclusive_variants(afm, t=.9)

            if one_dist:
                top_clone = n_vars.sort_values(ascending=False).index[0]
                n_variants_final[i] = n_vars.loc[top_clone]
                not_top_cells = labels.loc[lambda x: x != top_clone].index
                a_cells, _ = filter_cells_and_vars(
                    afm, cells=not_top_cells, variants=a_cells.var_names
                )
                final_labels.loc[labels.loc[lambda x: x == top_clone].index] = str(i)
                i += 1
            else:
                logger.info(f'Finished with {i} iterations: {t.stop()}. \n')

        except:
            logger.info('MQuad or vireo_vrapper ecountered errors. Finishing iterations.')
            break
    
    # Final labels
    logger.info(final_labels.value_counts())
    afm.obs = afm.obs.join(final_labels.to_frame('g'))

    ##

    # Process and viz
    df = afm.obs.join(
        final_labels.loc[lambda x: x.notna()].to_frame('MT_clones'),
        how='right'
    )
    df.to_csv(
        os.path.join(
            path_tmp, 
            f'take_one_off_{sample}_{filtering}_{n_muts}.csv'
        )
    )

    fig = contingency_iterative_plot(df, afm, good_quality_cells, figsize=(8,8))
    fig.tight_layout()
    fig.savefig(
        os.path.join(
            path_viz, f'take_one_off_{sample}_{filtering}_{n_muts}.png'
        ), 
        dpi=500
    )

    # Logs out
    logger.info(f'Finished job in {T.stop()}')

    ##


######################################################

# Run
if __name__ == '__main__':
    main()

