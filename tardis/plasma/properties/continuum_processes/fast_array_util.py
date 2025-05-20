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
        The result of cumulative integration of f along x
    """
    integ = (np.diff(x) * (f[1:] + f[:-1]) / 2.0).cumsum()
    return integ / integ[-1]


@njit(**njit_dict)
def cumulative_integrate_array_by_blocks(f, x, block_references):
    """
    Cumulatively integrate a function over blocks.

    This function cumulatively integrates a function `f` defined at
    locations `x` over blocks given in `block_references`.

    Parameters
    ----------
    f : numpy.ndarray, dtype float
        Input array to integrate. Shape is (N_freq, N_shells), where
        N_freq is the number of frequency values and N_shells is the number
        of computational shells.
    x : numpy.ndarray, dtype float
        The sample points corresponding to the `f` values. Shape is (N_freq,).
    block_references : numpy.ndarray, dtype int
        The start indices of the blocks to be integrated. Shape is (N_blocks,).

    Returns
    -------
    numpy.ndarray, dtype float
        Array with cumulatively integrated values. Shape is (N_freq, N_shells)
        same as f.
    """
    n_rows = len(block_references) - 1
    integrated = np.zeros_like(f)
    for i in prange(f.shape[1]):  # columns
        # TODO: Avoid this loop through vectorization of cumulative_trapezoid
        for j in prange(n_rows):  # rows
            start = block_references[j]
            stop = block_references[j + 1]
            integrated[start + 1 : stop, i] = numba_cumulative_trapezoid(
                f[start:stop, i], x[start:stop]
            )
    return integrated
