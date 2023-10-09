#!/usr/bin/env python
import numpy as np
from scipy.sparse import load_npz
from sklearn.preprocessing import normalize
import pickle

from fl_analysis.adaptive_walk import adaptive_walk

# parse Snakemake params
t_path = snakemake.input.transition_matrix
starting_set_path = snakemake.input.starting_genotypes
output = snakemake.output[0]
seed = int(snakemake.params.rseed)
sample_size = int(snakemake.params.sample_size)

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
data = adaptive_walk(T=T, local_peaks=local_peaks, sample_size=sample_size, max_steps=1000, starting_set=starting_set)

pickle.dump(data, open(output, "wb"))
