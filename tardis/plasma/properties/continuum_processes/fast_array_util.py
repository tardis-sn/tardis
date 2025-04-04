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
def numba_cumulative_trapezoid_2d(f, dx):
    """
    Cumulatively integrate f(x) using the composite trapezoidal rule.
    Vectorized trapezoidal integration for a single block.
    Precomputes `dx_block` externally to avoid redundant `np.diff`.

    Parameters
    ----------
    f : numpy.ndarray, dtype float
        Input array with multiple functions to integrate.
    dx : numpy.ndarray, dtype float
        The width of each block to integrate along.

    Returns
    -------
    numpy.ndarray, dtype float


    """
    n, m = f.shape
    trapz = dx.reshape(-1, 1) * (f[1:] + f[:-1]) / 2.0  # (n-1, m)
    integ = np.zeros_like(trapz)

    # Parallel cumulative sum over columns
    for i in prange(m):
        cumulative = 0.0
        for j in range(trapz.shape[0]):
            cumulative += trapz[j, i]
            integ[j, i] = cumulative
        if cumulative != 0:
            integ[:, i] /= cumulative
    return integ

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
    dx = np.diff(x)  # Precompute once globally
    integrated = np.zeros_like(f)
    n_blocks = len(block_references) - 1

    for j in prange(n_blocks):
        start = block_references[j]
        stop = block_references[j + 1]
        dx_block = dx[start: stop - 1]
        f_block = f[start:stop, :]
        integrated[start + 1:stop, :] = numba_cumulative_trapezoid_2d(f_block, dx_block)

    return integrated
