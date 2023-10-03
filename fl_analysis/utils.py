import numpy as np
import pandas as pd


def read_fl_data(path: str, delim: str):
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
    phe_arr = fl_df["score"].to_numpy(dtype=float)

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

