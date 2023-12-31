import numpy as np
import pandas as pd


# E. coli wild type transition & transversion biases according to Lee et al.
# PNAS, 109 (41), 2012 https://doi.org/10.1073/pnas.1210309109
mutation_biases = {('A', 'G'): .21,
                   ('T', 'C'): .21,
                   ('G', 'A'): .35,
                   ('C', 'T'): .35,
                   ('A', 'T'): .07,
                   ('T', 'A'): .07,
                   ('A', 'C'): .16,
                   ('T', 'G'): .16,
                   ('G', 'T'): .13,
                   ('C', 'A'): .13,
                   ('G', 'C'): .07,
                   ('C', 'G'): .07}


def add_mutation_bias(T, genotypes):
    for i, j in zip(*T.nonzero()):  # loop over sparse csr matrix
        g1, g2 = genotypes[i], genotypes[j]  # get genotypes
        mask = g1 != g2  # compare genotypes (find position where they differ)
        n1, n2 = g1[mask][0], g2[mask][0]  # get nucleotides
        T[(i, j)] *= mutation_biases[(n1, n2)]  # multiply fix. prob by mutation bias.

    return T


def get_column_from_csv(path: str, delim: str, col_id: str):
    """Read the data from the fitness landscape data file

    Args:
        path (str): Path to the fitness landscape data file
        delim (str): csv delimiter, e.g. "," or "\t"
        col_id (str): columns identifier
        
    Returns:
        np.array: column values as numpy array

    """
    fl_df = pd.read_csv(path, delimiter=delim)
    # turn into N x L array (N = #sequences, L = sequence length)
    arr = fl_df[col_id].to_numpy()
    return arr


def read_fl_data(path: str, delim: str, phe_column_id: str):
    """Read the data from the fitness landscape data file

    Args:
        path (str): Path to the fitness landscape data file
        delim (str): csv delimiter, e.g. "," or "\t"
        
    Returns:
        (np.array, np.array): array of genotypes, array of phenotypes 
    """
    fl_df = pd.read_csv(path, delimiter=delim)
    # turn into N x L array (N = #sequences, L = sequence length)
    seq_arr = np.array([list(seq) for seq in fl_df["seq"]])
    # remove columns that contain the same base, i.e. redundant sites of seq.
    seq_arr = remove_redundant(seq_arr)
    phe_arr = fl_df[phe_column_id].to_numpy(dtype=float)

    return seq_arr, phe_arr


def remove_redundant(arr):
    """Remove redundant columns from sequence array, i.e. an array with 
    dimensions NxL, where N is the number of sequences and L the seq. length.

    Args:
        arr (np.array): NxL array, N: number of sequences, L:
                                   sequence length.
    Returns:
        arr: np.ndarray
            output array without redundant columns
    """
    nonredundant_cols = []
    for i in range(arr.shape[1]):
        unique_elements = np.unique(arr[:, i])
        if len(unique_elements) > 1:
            nonredundant_cols.append(i)
    arr = arr[:, nonredundant_cols]
    return arr


def sequence_to_integers(sequence):
    """Turn an ATGC sequence into integer sequence.
    e.g. [[T, T, A, G, C],  -->  [[2, 2, 1, 3, 4],
          [A, G, C, T, T]]        [1, 3, 4, 2, 2]]

    Args:
        sequence (str): ATGC based sequence

    Returns:
        np.array: Array where bases (str) are replaced by integers

    """
    sequence_ = np.empty(sequence.shape, dtype=int)
    sequence_[sequence == 'A'] = 1
    sequence_[sequence == 'U'] = 2
    sequence_[sequence == 'G'] = 3
    sequence_[sequence == 'C'] = 4

    return sequence_


def kimura_fixation(p1: float, p2: float, N: int):
    """Formula to compute kimuara's fixation probability

    Args:
        p1 (float): Phenotype 1
        p2 (float): Phenotype 2
        N (int): Population size

    Returns:
        float: fixation probability (float in [0, 1]).
    
    """
    s = p2 - p1
    return (1-np.exp(-2*s))/(1-np.exp(-2*N*s))

