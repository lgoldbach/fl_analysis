#!/usr/bin/env python
import numpy as np

from fl_analysis.utils import n_highest_phenotypes

# parse Snakemake params
phenotype_path = snakemake.input.phenotypes
output = snakemake.output[0]
seed = int(snakemake.params.rseed)
starting_set_size = int(snakemake.params.starting_set_size)

# Set random seed
np.random.seed(seed)

ph = np.load(phenotype_path)

# Pick the n highest phenotypes
ph_high = n_highest_phenotypes(ph, n=starting_set_size)  

# ## this code is for picking random nodes that are not local peaks
# row_sums = np.array(T.sum(axis=1))[:,0]
# local_peaks = np.ma.masked_where(row_sums == 0, row_sums).mask

# non_peak_mask = np.invert(local_peaks)
# non_peaks = np.arange(T.shape[0])[non_peak_mask]

# starting_set = np.random.choice(non_peaks, size=starting_set_size, replace=False)  # pick 10000 random non-peak genotypes


np.save(output, ph_high)
