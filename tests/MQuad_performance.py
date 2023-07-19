#!/usr/bin/python

# MQuad timing test

# Code
import pickle
import sys
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


##


# Args
path_main = sys.argv[1]
sample = sys.argv[2]

# Paths
path_data = os.path.join(path_main, 'data') 
path_output = os.path.join(path_main, 'results', 'unsupervised_clones', 'output')
path_viz = os.path.join(path_main, 'results', 'unsupervised_clones', 'visualization')
path_tmp = os.path.join(path_main, 'results', 'unsupervised_clones', 'downstream_files')

# Set logger
path_ = os.getcwd()
make_folder(path_, 'logs', overwrite=False)
path_ = os.path.join(path_, 'logs')
logger = set_logger(path_, f'logs_MQuad_test.txt')


##


########################################################################

# Main
def main():

    t = Timer()

    ##

    # AFM first step, no filtering of low numebr of cells clones + variants filtering
    afm = read_one_sample(path_data, sample, with_GBC=True)

    # Here we go
    variants = {}
    for n in [100, 500, 1000, None]:

        # Filter and cluster
        t.start()
        logger.info(f'Beginning with {n} n_muts...')

        a_cells, a = filter_cells_and_vars(
            afm,
            filtering='MQuad', 
            min_cov_treshold=50,
            nproc=10,
            path_=os.getcwd(),
            n=n
        )
        variants[n] = a.var_names

        logger.info(f'Finished with {n} n_muts: {t.stop()}.')

    # Write down variants
    with open(os.path.join(path_tmp, 'test_MQuad.pickle'), 'wb') as f:
        pickle.dump(variants, f)


    ##


#########################################################################

# Run
if __name__ == '__main__':
    main()