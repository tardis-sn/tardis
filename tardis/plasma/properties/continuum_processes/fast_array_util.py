# It is currently not possible to use scipy.integrate.cumulative_trapezoid in
# numba. So here is my own implementation.
import numpy as np
from numba import njit, prange

from tardis.transport.montecarlo import njit_dict


@njit(**njit_dict)
def numba_cumulative_trapezoid(f, x):
    """
    Cumulatively integrate f(x) using the composite trapezoidal rule.

    Parameters
    ----------
    f : numpy.ndarray, dtype float
        Input array to integrate.
    x : numpy.ndarray, dtype float
        The coordinate to integrate along.

    Returns
    -------
    numpy.ndarray, dtype float
        The result of cumulative integration of f along x.
    """
    integ = (np.diff(x) * (f[1:] + f[:-1]) / 2.0).cumsum()
    result = np.zeros_like(f)  # Ensure the same size as `f`
    result[1:] = integ / integ[-1]  # Normalize and assign
    return result


@njit(**njit_dict)
def cumulative_integrate_array_by_blocks(f, x, block_references):
    """
    Perform cumulative integration over blocks using the trapezoidal rule, optimized.

    Parameters
    ----------
    f : numpy.ndarray, dtype float
        Input array to integrate (N_freq, N_shells).
    x : numpy.ndarray, dtype float
        Sample points corresponding to `f` (1D array, shared across blocks).
    block_references : numpy.ndarray, dtype int
        Start indices for each block (1D array, length N_blocks).

    Returns
    -------
    numpy.ndarray, dtype float
        Cumulatively integrated array (N_freq, N_shells).
    """
    integrated = np.zeros_like(f)
    n_blocks = len(block_references) - 1

    for block_idx in range(n_blocks):
        start = block_references[block_idx]
        stop = block_references[block_idx + 1]

        if stop - start > 1:  # Ensure the block has at least two points
            for col in range(f.shape[1]):  # Process each column (shell)
                integrated[start:stop, col] = numba_cumulative_trapezoid(
                    f[start:stop, col], x[start:stop]
                )

    return integrated
