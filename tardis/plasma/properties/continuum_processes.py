import logging

import numpy as np
import pandas as pd

from numba import prange, njit
from astropy import constants as const

from tardis.plasma.properties.base import ProcessingPlasmaProperty

__all__ = ["SpontRecombRateCoeff"]

logger = logging.getLogger(__name__)

njit_dict = {"fastmath": False, "parallel": False}


@njit(**njit_dict)
def integrate_array_by_blocks(f, x, block_references):
    """
    Integrates a function f defined at locations x over blocks
    given in block_references.

    Parameters
    ----------
    f : Two-dimensional Numpy Array, dtype float
    x : One-dimensional Numpy Array, dtype float
    block_references : One-dimensional Numpy Array, dtype int

    Returns
    -------
    integrated : Two-dimensional Numpy Array, dtype float

    """
    integrated = np.zeros((len(block_references) - 1, f.shape[1]))
    for i in prange(f.shape[1]):  # columns
        for j in prange(len(integrated)):  # rows
            start = block_references[j]
            stop = block_references[j + 1]
            integrated[j, i] = np.trapz(f[start:stop, i], x[start:stop])
    return integrated


def get_ion_multi_index(multi_index_full, next_higher=True):
    """
    Integrates a function f defined at locations x over blocks
    given in block_references.

    Parameters
    ----------
    multi_index_full : Pandas MultiIndex (atomic_number, ion_number,
                                          level_number)
    next_higher : bool
        If true use ion number of next higher ion, else use ion_number from
        multi_index_full.

    Returns
    -------
    multi_index : Pandas MultiIndex (atomic_number, ion_number)

    """
    atomic_number = multi_index_full.get_level_values(0)
    ion_number = multi_index_full.get_level_values(1)
    if next_higher is True:
        ion_number += 1
    return pd.MultiIndex.from_arrays([atomic_number, ion_number])


class SpontRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    alpha_sp : Pandas DataFrame, dtype float
               The rate coefficient for spontaneous recombination.
    """

    outputs = ("alpha_sp",)
    latex_name = (r"\alpha^{\textrm{sp}}",)

    def calculate(
        self,
        photo_ion_cross_sections,
        t_electrons,
        photo_ion_block_references,
        photo_ion_index,
        phi_ik,
    ):
        x_sect = photo_ion_cross_sections["x_sect"].values
        nu = photo_ion_cross_sections["nu"].values

        alpha_sp = 8 * np.pi * x_sect * nu ** 2 / (const.c.cgs.value) ** 2
        alpha_sp = alpha_sp[:, np.newaxis]
        boltzmann_factor = np.exp(
            -nu[np.newaxis].T
            / t_electrons
            * (const.h.cgs.value / const.k_B.cgs.value)
        )
        alpha_sp = alpha_sp * boltzmann_factor
        alpha_sp = integrate_array_by_blocks(
            alpha_sp, nu, photo_ion_block_references
        )
        alpha_sp = pd.DataFrame(alpha_sp, index=photo_ion_index)
        return alpha_sp * phi_ik
