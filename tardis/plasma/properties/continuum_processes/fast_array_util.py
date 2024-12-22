import numpy as np
from numba import njit, prange
from tardis.transport.montecarlo import njit_dict


@njit(**njit_dict)
def numba_cumulative_trapezoid(f, x):
    """
    Cumulatively integrate f(x) using the composite trapezoidal rule.

    Parameters
    ----------
    f : numpy.ndarray
        Input array to integrate. Shape: (N,)
    x : numpy.ndarray
        The coordinate array. Shape: (N,)

    Returns
    -------
    numpy.ndarray
        Cumulative integral of f along x, normalized by the final value. Shape: (N,)
    """
    if len(f) != len(x):
        raise ValueError("Input arrays f and x must have the same length.")
    if len(f) < 2:
        raise ValueError("Input arrays must have at least two elements for integration.")

    # Compute the cumulative trapezoidal integral
    dx = np.diff(x)
    cumulative = (dx * (f[1:] + f[:-1]) / 2.0).cumsum()

    # Normalize by the final value
    return np.concatenate(([0], cumulative / cumulative[-1]))


@njit(**njit_dict, parallel=True)
def cumulative_integrate_array_by_blocks(f, x, block_references):
    """
    Cumulatively integrate a 2D array over blocks defined by block references.

    Parameters
    ----------
    f : numpy.ndarray
        2D input array to integrate. Shape: (N_freq, N_shells)
    x : numpy.ndarray
        The coordinate array. Shape: (N_freq,)
    block_references : numpy.ndarray
        Start indices of the blocks to be integrated. Shape: (N_blocks,)

    Returns
    -------
    numpy.ndarray
        2D array with cumulative integrals for each block. Shape: (N_freq, N_shells)
    """
    n_blocks = len(block_references) - 1
    integrated = np.zeros_like(f)

    for col in prange(f.shape[1]):  # Iterate over columns (N_shells)
        for block_idx in range(n_blocks):  # Iterate over blocks
            start = block_references[block_idx]
            stop = block_references[block_idx + 1]

            if stop - start < 2:
                continue  # Skip blocks that are too small to integrate

            integrated[start:stop, col] = numba_cumulative_trapezoid(
                f[start:stop, col], x[start:stop]
            )

    return integrated
