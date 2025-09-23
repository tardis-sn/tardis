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
<<<<<<< HEAD
@njit(**njit_dict) 
=======


@njit(**njit_dict)
>>>>>>> 46e2fe9f1c6f954b93947c9fbc4f638f3e261f5e
def cumulative_integrate_array_by_blocks(f, x, block_references, normalize=True):
    """
    Cumulatively integrates an array f(x) over discrete blocks.

    This function performs a trapezoidal integration. The calculation is done
    independently for each block defined in block_references.

    Parameters
    ----------
    f : np.ndarray
        The function values, with shape (N, M).
    x : np.ndarray
        The coordinates, with shape (N,).
    block_references : np.ndarray
        An array of indices defining the start and end of each block.
        e.g., [0, 10, 25] defines two blocks: [0:10] and [10:25].
    normalize : bool, optional
        If True, the cumulative integral within each block is normalized
        so that its final value is 1. Defaults to True.

    Returns
    -------
    np.ndarray
        An array with the same shape as f containing the cumulatively
        integrated values within each block.
    """
    integrated = np.zeros_like(f)
<<<<<<< HEAD
    n_blocks = len(block_references) - 1

    if n_blocks <= 0:
        return integrated 

    for i in prange(f.shape[1]):
       
        for j in prange(n_blocks):
            start = block_references[j]
            stop = block_references[j + 1]
            block_len = stop - start

            if block_len < 2:
                if block_len == 1:
                    integrated[start, i] = f[start, i]
                continue

            f_block = f[start:stop, i]
            x_block = x[start:stop]

            dx = np.diff(x_block)
            contribs = dx * 0.5 * (f_block[1:] + f_block[:-1])

            block_integral = np.zeros(block_len)
            block_integral[1:] = np.cumsum(contribs)

            if normalize and block_integral[-1] != 0:
                block_integral /= block_integral[-1]
            integrated[start:stop, i] = block_integral
            
=======
    dx = np.diff(x)

    for i in prange(f.shape[1]):
        contribs = dx * 0.5 * (f[1:, i] + f[:-1, i])
        cumsum = np.zeros(f.shape[0])
        cumsum[1:] = np.cumsum(contribs)

        for j in range(len(block_references) - 1):
            start = block_references[j]
            stop = block_references[j + 1]

            block_len = stop - start
            if block_len == 1:
                # Single-element block: copy original value
                integrated[start, i] = f[start, i]
            else:
                block = cumsum[start:stop] - cumsum[start]
                if normalize and block[-1] != 0:
                    block /= block[-1]
                # Fill the integrated array, aligned properly
                integrated[start + 1:stop, i] = block[1:]  # skip first 0

>>>>>>> 46e2fe9f1c6f954b93947c9fbc4f638f3e261f5e
    return integrated