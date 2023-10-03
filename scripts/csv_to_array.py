#!/usr/bin/env python
import numpy as np

from fl_analysis.utils import read_fl_data, remove_redundant, sequence_to_integers

# Snakemake input output definitions
input = snakemake.input[0]
output_genotypes = snakemake.output.genotypes
output_genotypes_num = snakemake.output.genotypes_num
output_phenotypes = snakemake.output.phenotypes
csv_delimiter = snakemake.params.csv_delimiter

seq_arr, phe_arr = read_fl_data(input, delim=csv_delimiter)
seq_arr_num = sequence_to_integers(seq_arr)  # encode bases as integers

np.save(output_genotypes, seq_arr)
np.save(output_genotypes_num, seq_arr_num)
np.save(output_phenotypes, phe_arr)
