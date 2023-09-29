#!/usr/bin/env python
import numpy as np
from scipy.sparse import load_npz
from sklearn.preprocessing import normalize
import pickle

from fl_analysis.adaptive_walk import adaptive_walk

# parse Snakemake params
t_path = snakemake.input[0]
output = snakemake.output[0]
seed = int(snakemake.params.rseed)
starting_set_size = int(snakemake.params.starting_set_size)
sample_size = int(snakemake.params.sample_size)

# Set random seed
np.random.seed(seed)

T = load_npz(t_path)
Tn = T.copy()
Tn = normalize(Tn, norm="l1")

row_sums = np.array(Tn.sum(axis=1))[:,0]
local_peaks = np.ma.masked_where(row_sums == 0, row_sums).mask

non_peak_mask = np.invert(local_peaks)
non_peaks = np.arange(Tn.shape[0])[non_peak_mask]
starting_set_size = len(non_peaks)
starting_set = np.random.choice(non_peaks, size=starting_set_size, replace=False)  # pick 10000 random non-peak genotypes

data = adaptive_walk(T=Tn, local_peaks=local_peaks, sample_size=sample_size, max_steps=1000, starting_set=starting_set)
pickle.dump(data, open(output, "wb"))
