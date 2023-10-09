#!/usr/bin/env python
import numpy as np
from scipy.sparse import csr_matrix, save_npz
from sklearn.preprocessing import normalize

from fl_analysis.utils import kimura_fixation, add_mutation_bias


# Snakemake input output definitions
input_genotypes = snakemake.input.genotypes
input_phenotypes = snakemake.input.phenotypes
input_hamming_distances = snakemake.input.hamming_dist
output = snakemake.output[0]
population_size = 10 ** int(snakemake.params.pop_size_exp)

hamming_dist = np.load(input_hamming_distances)
ph1 = np.load(input_phenotypes[0])
ph2 = np.load(input_phenotypes[1])

edges = np.where(hamming_dist == 1)  # get edges of neighbors

# compute differences between phenotypes for both phenotypes in both directions
ph_diff1_1 = ph1[edges[0]] - ph1[edges[1]]
ph_diff1_2 = ph1[edges[1]] - ph1[edges[0]]
ph_diff2_1 = ph2[edges[0]] - ph2[edges[1]]
ph_diff2_2 = ph2[edges[1]] - ph2[edges[0]]

# overall phenotype difference between two genotypes is the sum of both 
# phenotype differences
s1 = ph_diff1_1 + ph_diff2_1  
s2 = ph_diff1_2 + ph_diff2_2

kimura_fixation_ = lambda s: kimura_fixation(s, N=population_size)
kimura_vect = np.vectorize(kimura_fixation_)  # vectorize transition prob. func.

# compute fixation probability aka. transition prob. in both directions (1, 2)
fix_prob1 = kimura_vect(s1)
fix_prob1 = np.nan_to_num(fix_prob1)
fix_prob2 = kimura_vect(s2)
fix_prob2 = np.nan_to_num(fix_prob2)

fix_probs = np.concatenate((fix_prob1, fix_prob2))

edges0 = np.concatenate((edges[1], edges[0]))
edges1 = np.concatenate((edges[0], edges[1]))

dim = ph1.shape[0]
T = csr_matrix((fix_probs, (edges0, edges1)), shape=(dim, dim))
T = add_mutation_bias(T=T, genotypes=np.load(input_genotypes))
T = normalize(T, norm="l1")

save_npz(output, T)
