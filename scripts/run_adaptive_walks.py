#!/usr/bin/env python
import numpy as np
from scipy.sparse import load_npz
from sklearn.preprocessing import normalize
import pickle
from snakemake.script import Snakemake
import argparse

from fl_analysis.adaptive_walk import adaptive_walk

try:
    # check if script is called from snakemake
    if isinstance(snakemake, Snakemake):  
        # parse Snakemake params
        t_path = snakemake.input.transition_matrix
        starting_set_path = snakemake.input.starting_genotypes
        output = snakemake.output[0]
        seed = int(snakemake.params.rseed)
        max_steps = int(snakemake.params.max_steps)
        sample_size = int(snakemake.params.sample_size)
        
except NameError:
    # or from command-line
    if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("-t", "--tmatrix", help="Transition matrix", type=str)
        parser.add_argument("-s", "--starting_set", help="Set of starting genotypes", type=str)
        parser.add_argument("-o", "--output", help="Output file", type=str)
        parser.add_argument("-r", "--seed", help="Random seed", type=int)
        parser.add_argument("-m", "--max_steps", help="Maximum number of steps", type=int)
        parser.add_argument("-n", "--sample_size", help="Number of paths to sample per starting genotype", type=int)

        args = parser.parse_args()

        t_path = args.tmatrix
        starting_set_path = args.starting_set
        output = args.output
        seed = args.seed
        max_steps = args.max_steps
        sample_size = args.sample_size
        

# Set random seed
np.random.seed(seed)

T = load_npz(t_path)
starting_set = np.load(starting_set_path)

row_sums = np.array(T.sum(axis=1))[:,0]
# local_peaks = np.ma.masked_where(row_sums == 0, row_sums, shrink=False).mask

# get mask that is True where local peaks are (cannot use np.ma.masked_where 
# because that returns only False on the case where we have no local peaks, 
# e.g low population sizes)
ma = np.ma.masked_equal(row_sums, 0)
local_peaks = np.ma.getmaskarray(ma)
data = adaptive_walk(T=T, local_peaks=local_peaks, sample_size=sample_size, max_steps=max_steps, starting_set=starting_set)

pickle.dump(data, open(output, "wb"))
