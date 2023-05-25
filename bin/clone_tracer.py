#!/usr/bin/python

# Clonal inference with CloneTracer

########################################################################

# Code
import argparse
import os

##

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='clone_tracer',
    description='Wrapper script within CloneTracer run_clonetracer.py script.'
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

# Probability treshold
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

# GPU
my_parser.add_argument(
    '-g', 
    '--with_gpu', 
    default=False, 
    action='store_true',
    help='Boolean indicating whether to run w/i or w/o GPU. Default: False'
)

# Parse arguments
args = my_parser.parse_args()


##


path_data = args.path_data
path_ = os.getcwd()
sample = args.sample
filtering = args.filtering
min_cov_treshold = args.min_cov_treshold
min_cell_number = args.min_cell_number
p_treshold = args.p_treshold
ncores = args.ncores
gpu = args.with_gpu


path_data = '/Users/IEO5505/Desktop/mito_bench/data'
sample = 'AML_clones'
path_ = os.getcwd()
filtering = 'miller2022'
min_cov_treshold = 50
min_cell_number = 10
ncores = 8
gpu = False


##


########################################################################

# Code
import pickle
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
from mito_utils.embeddings_plots import *
from mito_utils.dimred import *
from helper_functions import *
import warnings
warnings.filterwarnings("ignore")
pyro.set_rng_seed(1234)

####################################################################

# Main
def main():

    # Load data
    T = Timer()
    T.start()
    
    t = Timer()

    # Logging
    path_ = os.getcwd()
    logger = set_logger(path_, f'log_CloneTracer_{sample}_{filtering}_{min_cell_number}.txt')
    
    # Set logger
    logger.info(
        f""" 
        Execute CloneTracer: \n
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

    # Filter cells and vars, create a mutational embedding
    _, a = filter_cells_and_vars(
        afm,
        sample=sample,
        filtering=filtering, 
        min_cell_number=min_cell_number,
        min_cov_treshold=min_cov_treshold,
        nproc=ncores, 
        path_=path_
    )

    # Extract filtered feature matrix, format and reduce with UMAP
    a = nans_as_zeros(a) # For sklearn APIs compatibility

    # Get parallel matrices
    AD, DP, _ = get_AD_DP(a, to='csc')
    
    # Prep input 
    d_input = {
        'M' : AD.A.T.tolist(), # Transpose here
        'N' : DP.A.T.tolist(),
        'mut_type' : [ 2 for _ in a.var_names ],                   # Only MT-SNPs here
        'mut_names' : list(a.var_names),
        'cell_barcode' : list(a.obs_names)
    }
    logger.info(f'Input preparation {t.stop()}')

    # Set tree instance
    torch.set_default_tensor_type(torch.DoubleTensor)
    
    
    
import json
with open(input_file) as f:
    d = json.load(f)
 
    

    input_file = '/Users/IEO5505/CloneTracer/clonal_inference/data/A.6.json'
    name = 'A.6'
    mult_samp = True
    cnv_celltype = False
    gpu = False

    t = create_tree_class(input_file, 'A.6', mult_samp=True, cnv_celltype=False, gpu=False)
    

    

    # Run SVI
    logger.info(f'Begin inference...')
    t.infer_hierarchy(200, 190, os.getcwd()) # Da sistemare num_iter e init
    logger.info(f'Finished inference: {t.stop()}')
    
    # Write output
    with open(f'CloneTracer_{sample}.pickle', 'wb') as f:
       pickle.dump(t, f)
       
    # Exit
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

########################################################################

# Run program
if __name__ == "__main__":
    main()
    
    
    
    
    
    
# function to create tree class object from JSON input file
def create_tree_class(input_file, name, mult_samp, cnv_celltype, gpu):
    

    # load data from JSON file
    with open(input_file, "rb") as f:
        input_data = json.load(f)

    # get number of mutations
    nmuts = len(input_data["mut_names"])

    # if there are more than one sample

    # make a dictionary with input data
    data_svi = {"M": torch.Tensor(input_data["M"]),
                 "N": torch.Tensor(input_data["N"]),
                 "mut_type": torch.Tensor(input_data["mut_type"]),
                 "names": input_data["mut_names"],
                 "barcodes": input_data["cell_barcode"],
                 "class_af": mult_samp,
                 "cnv_celltype": cnv_celltype}

    # bulk data if present
    for entry in ["bulk_M", "bulk_N", "r_cnv"]:

        # if present add information to dictionary
        if entry in input_data and input_data[entry]:

            data_svi[entry] = torch.Tensor(input_data[entry])
            
            if entry == "bulk_M":
                
                data_svi["bulk_af"] = True

        # otherwise set values to 0
        else:
            data_svi[entry] = torch.zeros(nmuts)
            
            if entry == "bulk_M":
                
                data_svi["bulk_af"] = False

    # priors for heteroplasmy (default values are 1000,1000 for nuclear, 1,1 for mitochondria and 2,100 for CNVs)
    for entry in ["h_alpha", "h_beta"]:

        # if present add information to dictionary
        if entry in input_data and input_data[entry]:

            data_svi[entry] = torch.Tensor(input_data[entry])

        # otherwise set to default values
        else:

            if entry == "h_alpha":

                h_mapper = {0: 2, 1: 1000, 2: 1}

                data_svi[entry] = torch.Tensor([h_mapper[mut] for mut in input_data["mut_type"]])

            else:

                h_mapper = {0: 100, 1: 1000, 2: 1}

                data_svi[entry] = torch.Tensor([h_mapper[mut] for mut in input_data["mut_type"]])
        
    # convert dictionary with cnv priors to tensor    
    if cnv_celltype:
        
        mean = list()
        sd = list()
        
        for chrom in input_data["cnv_priors"]:
            
            mean.append(input_data["cnv_priors"][chrom]["mean"])
            sd.append(input_data["cnv_priors"][chrom]["sd"])
            
        data_svi["cnv_ct_mean"] = torch.Tensor(mean)
        data_svi["cnv_ct_sd"] = torch.Tensor(sd)
        
    else:
        
        data_svi["cnv_ct_mean"] = []
        data_svi["cnv_ct_sd"] = []

    # add additional information for celltype-specific CNV model (if present)
    for entry in ["class_assign", "class_names", "celltype", "celltype_names", "umapx", "umapy"]:

        if entry in input_data and input_data[entry]:

            if entry in ["class_assign"] and gpu:

                data_svi[entry] = torch.cuda.IntTensor(input_data[entry])

            elif entry in ["class_assign"] and not gpu:

                data_svi[entry] = torch.IntTensor(input_data[entry])
                
            elif entry in ["celltype"]:
                
                data_svi[entry] = torch.tensor(input_data[entry])
                
            else:

                data_svi[entry] = input_data[entry]

        else:

            data_svi[entry] = []

    # rename bulk data entries
    data_svi["af_alpha"] = data_svi["bulk_M"]
    data_svi["af_beta"] = data_svi["bulk_N"]

    t = tree(name, data_svi)

    return(t)