#!/usr/bin/env python
import numpy as np
from scipy.spatial.distance import hamming

# Snakemake input output definitions
input = snakemake.input[0]
output = snakemake.output[0]

# read data and create empty array for hamming distances
seq_arr_num = np.load(input)
num_of_seq = seq_arr_num.shape[0]
seq_len = seq_arr_num.shape[1]
hamming_dist = np.empty((num_of_seq, num_of_seq))

# compute pairwise hamming distances
for i in range(num_of_seq):
    for j in range(i):
        
        hamming_dist[i, j] = hamming(seq_arr_num[i], seq_arr_num[j]) 

# multiply by sequence length to get hamming distance not as fraction
# but as mutational step. Default of scipy hamming is fraction [0, 1].
hamming_dist * seq_len

np.save(output, hamming_dist)
