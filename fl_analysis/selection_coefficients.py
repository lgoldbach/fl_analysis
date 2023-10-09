import numpy as np
import pandas as pd


def phenotype_diff(ph_current: np.array, ph_mut: np.array):
    """Compute the selection coefficient between two arrays of phenotypes as 
    the difference between the two

    Args:
        ph_current (1D np.array): Phenotype array before mutation
        ph_mut (1D np.array): Phenotype array after mutation
    
    Returns:
        np.array: array of selection coefficients

    """
    return ph_current - ph_mut


def double_phenotype_diff(ph1_current, ph1_mut, ph2_current, ph2_mut):
    """Compute the selection coefficient based on two different phenotypes

    Args:
        ph1_current (1D np.array): Phenotype array before mutation
        ph1_current (1D np.array): Phenotype array after mutation
        ph2_current (1D np.array): Phenotype array before mutation
        ph2_current (1D np.array): Phenotype array after mutation
    
    Returns:
        np.array: array of selection coefficients

    """
    s1 = ph1_current - ph1_mut
    s2 = ph2_current - ph1_mut
    return s1 + s2
