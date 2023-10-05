#!/usr/bin/env python
import numpy as np

from fl_analysis.utils import get_column_from_csv, remove_redundant, sequence_to_integers

# Snakemake input output definitions
input = snakemake.input[0]
output_phenotypes = snakemake.output.phenotype
csv_delimiter = snakemake.params.csv_delimiter
csv_pheno_id = snakemake.params.csv_pheno_id

ph = get_column_from_csv(input, delim=csv_delimiter, col_id=csv_pheno_id)

np.save(output_phenotypes, ph)
