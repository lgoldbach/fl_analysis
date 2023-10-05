#!/usr/bin/env python
import numpy as np

from fl_analysis.utils import get_column_from_csv, remove_redundant, sequence_to_integers

# Snakemake input output definitions
input = snakemake.input[0]
output_genotypes = snakemake.output.genotypes
output_genotypes_num = snakemake.output.genotypes_num
csv_delimiter = snakemake.params.csv_delimiter

seq_arr = get_column_from_csv(input, delim=csv_delimiter, col_id="seq")
# turn genotypes strings to list of letters
seq_arr = np.array([list(seq) for seq in seq_arr])
# remove columns that contain the same base, i.e. redundant sites of seq.
seq_arr = remove_redundant(seq_arr) 
seq_arr_num = sequence_to_integers(seq_arr)  # encode bases as integers

np.save(output_genotypes, seq_arr)
np.save(output_genotypes_num, seq_arr_num)
