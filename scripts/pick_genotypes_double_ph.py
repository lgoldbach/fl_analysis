#!/usr/bin/env python
import numpy as np
from fl_analysis.utils import n_lowest_phenotypes

# parse Snakemake params
genotypes = snakemake.input.genotypes
output = snakemake.output[0]
seed = int(snakemake.params.rseed)
starting_set_size = int(snakemake.params.starting_set_size)

# Get runs for all genotypes
gt = np.load(genotypes)
gt_ids = np.arange(gt.shape[0])  

np.save(output, gt_ids)


# # Set random seed
# np.random.seed(seed)

# ph2 = np.load(phenotype_paths[1])

# # pick the phenotype that has the lowest combined distance from the 
# # respective mean phenotype
# ph1 = np.abs(ph1 - np.mean(ph1))
# ph2 = np.abs(ph2 - np.mean(ph2))

# ph = ph1 + ph2

# ph_ids = n_lowest_phenotypes(ph, n=starting_set_size)  


