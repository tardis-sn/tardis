"""Import-light kernels for homologous Sobolev optical depths."""

import numpy as np
import numpy.typing as npt
import pandas as pd
from astropy import units as u
from numba import jit, prange

from tardis import constants as const

SOBOLEV_COEFFICIENT = (
    (
        ((np.pi * const.e.gauss**2) / (const.m_e.cgs * const.c.cgs))
        * u.cm
        * u.s
        / u.cm**3
    )
    .to(1)
    .value
)


def calculate_sobolev_line_opacity(
    lines: pd.DataFrame,
    level_number_density: pd.DataFrame,
    time_explosion: u.Quantity,
    stimulated_emission_factor: npt.NDArray[np.float64],
) -> pd.DataFrame:
    """Calculate homologous Sobolev optical depths for all lines.

    Parameters
    ----------
    lines : pandas.DataFrame
        Atomic line data.
    level_number_density : pandas.DataFrame
        Level populations indexed by atomic level and columned by shell.
    time_explosion : astropy.units.Quantity
        Time since explosion.
    stimulated_emission_factor : numpy.ndarray
        Stimulated-emission correction indexed by line and shell.

    Returns
    -------
    pandas.DataFrame
        Sobolev optical depth indexed by line and columned by shell.
    """
    tau_sobolevs = (
        (lines.wavelength_cm * lines.f_lu).values[np.newaxis].T
        * SOBOLEV_COEFFICIENT
        * time_explosion.to(u.s).value
        * stimulated_emission_factor
        * level_number_density.reindex(lines.droplevel(-1).index).values
    )
    if np.any(np.isnan(tau_sobolevs)) or np.any(np.isinf(np.abs(tau_sobolevs))):
        raise ValueError(
            "Some tau_sobolevs are nan, inf, -inf in tau_sobolevs."
            " Something went wrong!"
        )
    return pd.DataFrame(
        tau_sobolevs,
        index=lines.index,
        columns=np.array(level_number_density.columns),
    )


@jit(nopython=True, parallel=True)
def numba_calculate_beta_sobolev(
    tau_sobolevs: npt.NDArray[np.float64],
    beta_sobolevs: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Calculate escape probabilities in place for flattened optical depths."""
    for i in prange(len(tau_sobolevs)):
        if tau_sobolevs[i] > 1e3:
            beta_sobolevs[i] = tau_sobolevs[i] ** -1
        elif tau_sobolevs[i] < 1e-4:
            beta_sobolevs[i] = 1 - 0.5 * tau_sobolevs[i]
        else:
            beta_sobolevs[i] = (1 - np.exp(-tau_sobolevs[i])) / (
                tau_sobolevs[i]
            )
    return beta_sobolevs


def calculate_beta_sobolev(tau_sobolevs: pd.DataFrame) -> pd.DataFrame:
    """Calculate Sobolev escape probabilities from optical depths.

    Parameters
    ----------
    tau_sobolevs : pandas.DataFrame
        Sobolev optical depths indexed by line and columned by shell.

    Returns
    -------
    pandas.DataFrame
        Sobolev escape probabilities with the same labels and shape.
    """
    tau_values_1d = tau_sobolevs.to_numpy().ravel()
    beta_values_1d = np.zeros(tau_values_1d.shape[0], dtype=np.float64)
    numba_calculate_beta_sobolev(tau_values_1d, beta_values_1d)
    return pd.DataFrame(
        beta_values_1d.reshape(tau_sobolevs.shape),
        index=tau_sobolevs.index,
        columns=tau_sobolevs.columns,
    )
