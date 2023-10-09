#!/usr/bin/env python
import numpy as np
from scipy.sparse import csr_matrix, save_npz
from sklearn.preprocessing import normalize

from fl_analysis.utils import add_mutation_bias

# Snakemake input output definitions
input_genotypes = snakemake.input.genotypes
input_phenotypes = snakemake.input.phenotypes
input_hamming_distances = snakemake.input.hamming_dist
output = snakemake.output[0]

hamming_dist = np.load(input_hamming_distances)
ph = np.load(input_phenotypes)
edges = np.where(hamming_dist == 1)  # get edges of neighbors

# compute differences between phenotypes
ph_diff1 = ph[edges[0]] - ph[edges[1]]
ph_diff2 = ph[edges[1]] - ph[edges[0]]
ph_diff1[ph_diff1 < 0] = 0  # set all differences below 0 to 0 because we do not care about fitness-descending steps    
ph_diff2[ph_diff2 < 0] = 0

score_diffs = np.concatenate((ph_diff1, ph_diff2))

edges0 = np.concatenate((edges[1], edges[0]))  # get both directions for each edge
edges1 = np.concatenate((edges[0], edges[1]))

dim = ph.shape[0]
T = csr_matrix((score_diffs, (edges0, edges1)), shape=(dim, dim))
T.eliminate_zeros()
T[T.nonzero()] = 1
T = add_mutation_bias(T=T, genotypes=np.load(input_genotypes))
T = normalize(T, norm="l1")

save_npz(snakemake.output[0], T)
