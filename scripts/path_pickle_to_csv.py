#!/usr/bin/env python
import numpy as np
import pickle
import pandas as pd


# parse Snakemake params
genotypes_path = snakemake.input.genotypes
paths_path = snakemake.input.paths
output = snakemake.output[0]

paths_raw, counts = pickle.load(open(paths_path, "rb"))
genotypes = np.load(genotypes_path)

# turn into dict
paths = {}
for p in paths_raw:
    paths[p] = counts[paths_raw[p]]

df1 = pd.DataFrame(paths.items(), columns=["path", "count"])
df1.to_csv(output, index=False)
