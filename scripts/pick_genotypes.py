#!/usr/bin/env python
import numpy as np

from fl_analysis.utils import n_highest_phenotypes, n_lowest_phenotypes

# parse Snakemake params
phenotype_path = snakemake.input.phenotypes
output = snakemake.output[0]
seed = int(snakemake.params.rseed)
starting_set_size = int(snakemake.params.starting_set_size)
# by which criterion to pick genotypes 
# ("lowest" - n lowest phenotypes, "highest" - n highest phenotypes, "random")
criterion = snakemake.params.criterion

# Set random seed
np.random.seed(seed)

ph = np.load(phenotype_path)

if criterion == "highest":
    # Pick the n highest phenotypes
    ph_ids = n_highest_phenotypes(ph, n=starting_set_size)  
elif criterion == "lowest":
    # Pick the n highest phenotypes
    ph_ids = n_lowest_phenotypes(ph, n=starting_set_size)
elif criterion == "random":
    ph_ids = np.random.choice(ph.size, size=starting_set_size)

    # ## this code is for picking random nodes that are not local peaks
    # T = np.load(transition_matrix)
    # row_sums = np.array(T.sum(axis=1))[:,0]
    # local_peaks = np.ma.masked_where(row_sums == 0, row_sums).mask

    # non_peak_mask = np.invert(local_peaks)
    # non_peaks = np.arange(T.shape[0])[non_peak_mask]

    # starting_set = np.random.choice(non_peaks, size=starting_set_size, replace=False)  # pick 10000 random non-peak genotypes

np.save(output, ph_ids)
