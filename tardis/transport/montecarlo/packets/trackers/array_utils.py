"""
Array utility functions for packet tracking.

These functions provide common array extension operations used by
the packet tracking system.
"""
import numpy as np
from numba import njit

from tardis.transport.montecarlo.packets.radiative_packet import InteractionType

NO_INTERACTION_INT = int(InteractionType.NO_INTERACTION)


@njit
def extend_array(array: np.ndarray, new_length: int) -> np.ndarray:
    """
    Extend array to new length preserving existing data.

    Parameters
    ----------
    array : np.ndarray
        Array to extend.
    new_length : int
        New length for the array.

    Returns
    -------
    np.ndarray
        Extended array with original data preserved.
    """
    temp_array = np.empty(new_length, dtype=array.dtype)
    temp_array[:len(array)] = array
    return temp_array


@njit
def extend_interaction_type_array(array: np.ndarray, new_length: int) -> np.ndarray:
    """
    Extend interaction type array with NO_INTERACTION default.

    Parameters
    ----------
    array : np.ndarray
        Array to extend.
    new_length : int
        New length for the array.

    Returns
    -------
    np.ndarray
        Extended array with NO_INTERACTION defaults for new elements.
    """
    temp_array = np.full(new_length, NO_INTERACTION_INT, dtype=array.dtype)
    temp_array[:len(array)] = array
    return temp_array


@njit
def extend_float_array(array: np.ndarray, new_length: int) -> np.ndarray:
    """
    Extend float array with NaN default.

    Parameters
    ----------
    array : np.ndarray
        Array to extend.
    new_length : int
        New length for the array.

    Returns
    -------
    np.ndarray
        Extended array with NaN defaults for new elements.
    """
    temp_array = np.full(new_length, np.nan, dtype=array.dtype)
    temp_array[:len(array)] = array
    return temp_array


@njit
def extend_int_array(array: np.ndarray, new_length: int) -> np.ndarray:
    """
    Extend integer array with -1 default.

    Parameters
    ----------
    array : np.ndarray
        Array to extend.
    new_length : int
        New length for the array.

    Returns
    -------
    np.ndarray
        Extended array with -1 defaults for new elements.
    """
    temp_array = np.full(new_length, -1, dtype=array.dtype)
    temp_array[:len(array)] = array
    return temp_array
