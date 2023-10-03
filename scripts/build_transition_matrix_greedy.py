#!/usr/bin/env python
import numpy as np
from scipy.sparse import csr_matrix, save_npz
from sklearn.preprocessing import normalize


# Snakemake input output definitions
input_phenotypes = snakemake.input.phenotypes
input_hamming_distances = snakemake.input.hamming_dist
output = snakemake.output[0]

hamming_dist = np.load(input_hamming_distances)
ph = np.load(input_phenotypes)

edges = np.where(hamming_dist == 1)  # get edges of neighbors

# compute differences between phenotypes
ph_diff1 = ph[edges[0]] - ph[edges[1]]
ph_diff2 = ph[edges[1]] - ph[edges[0]]

score_diffs = np.concatenate((ph_diff1, ph_diff2))

edges0 = np.concatenate((edges[1], edges[0]))  # get both directions for each edge
edges1 = np.concatenate((edges[0], edges[1]))

dim = ph.shape[0]
T = csr_matrix((score_diffs, (edges0, edges1)), shape=(dim, dim))


# find the phenotype maximizing moves and construct a new matrix from these
# faster than changing csr_matrix in place
data = []
rows = []
cols = []
for i, row in enumerate(T):
    max_i = row.argmax()  # get index of max val
    max = T[i, max_i]  # get max. val
    if max > 0:  # if there is no higher neighbor
        rows.append(i)
        cols.append(max_i)
        data.append(1)  # set maximizing move to prob 1
T = csr_matrix((data, (rows, cols)), shape=(dim, dim))

save_npz(snakemake.output[0], T)
