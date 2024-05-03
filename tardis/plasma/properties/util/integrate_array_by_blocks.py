import numpy as np
from numba import njit, prange

from tardis.transport.montecarlo import njit_dict


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
