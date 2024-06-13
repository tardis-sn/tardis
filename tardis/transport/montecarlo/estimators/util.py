import numpy as np
import pandas as pd
from numba import njit, prange

from tardis.transport.montecarlo import njit_dict


def bound_free_estimator_array2frame(
    bound_free_estimator_array, level2continuum_idx
):
    """
    Transform a bound-free estimator array to a DataFrame.

    This function transforms a bound-free estimator array with entries
    sorted by frequency to a multi-indexed DataFrame sorted by level.

    Parameters
    ----------
    bf_estimator_array : numpy.ndarray, dtype float
        Array of bound-free estimators (e.g., for the stimulated recombination rate)
        with entries sorted by the threshold frequency of the bound-free continuum.
    level2continuum_idx : pandas.Series, dtype int
        Maps a level MultiIndex (atomic_number, ion_number, level_number) to
        the continuum_idx of the corresponding bound-free continuum (which are
        sorted by decreasing frequency).

    Returns
    -------
    pandas.DataFrame, dtype float
        Bound-free estimators indexed by (atomic_number, ion_number, level_number).
    """
    bf_estimator_frame = pd.DataFrame(
        bound_free_estimator_array, index=level2continuum_idx.index
    ).sort_index()
    bf_estimator_frame.columns.name = "Shell No."
    return bf_estimator_frame


@njit(**njit_dict)
def integrate_array_by_blocks(f, x, block_references):
    """
    Integrate a function over blocks.

    This function integrates a function `f` defined at locations `x`
    over blocks given in `block_references`.

    Parameters
    ----------
    f : numpy.ndarray, dtype float
        2D input array to integrate.
    x : numpy.ndarray, dtype float
        1D array with the sample points corresponding to the `f` values.
    block_references : numpy.ndarray, dtype int
        1D array with the start indices of the blocks to be integrated.

    Returns
    -------
    numpy.ndarray, dtype float
        2D array with integrated values.
    """
    integrated = np.zeros((len(block_references) - 1, f.shape[1]))
    for i in prange(f.shape[1]):  # columns
        for j in prange(len(integrated)):  # rows
            start = block_references[j]
            stop = block_references[j + 1]
            integrated[j, i] = np.trapz(f[start:stop, i], x[start:stop])
    return integrated
