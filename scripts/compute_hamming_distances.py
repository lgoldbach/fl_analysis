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
hamming_dist = np.zeros((num_of_seq, num_of_seq), dtype=int)

# compute pairwise hamming distances
for i in range(num_of_seq):
    for j in range(i):
        # multiply by sequence length to get hamming distance not as fraction
        # but as mutational step. Default of scipy hamming is fraction [0, 1].
        hamming_dist[i, j] = int(hamming(seq_arr_num[i], seq_arr_num[j]) * seq_len)

np.save(output, hamming_dist)
