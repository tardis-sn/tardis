import numpy as np
import pandas as pd
from astropy import units as u
from numba import jit, prange

from tardis import constants as const
from tardis.plasma.properties.base import ProcessingPlasmaProperty

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
    lines,
    level_number_density,
    time_explosion,
    stimulated_emission_factor,
):
    """
    Calculates the Sobolev line opacity based on the provided parameters.

    Parameters
    ----------
    lines : pandas.DataFrame
        DataFrame containing information about spectral lines.
    level_number_density : pandas.DataFrame
        DataFrame with level number densities.
    time_explosion : astropy.units.Quantity
        Time since explosion.
    stimulated_emission_factor : float
        Factor for stimulated emission.

    Returns
    -------
    pandas.DataFrame
        Calculated Sobolev line opacity values.

    Raises
    ------
    ValueError
        If any calculated tau_sobolevs are nan or inf.

    Examples
    --------
    >>> calculate_sobolev_line_opacity(lines_data, level_density_data, time_exp, stim_factor)
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
def numba_calculate_beta_sobolev(tau_sobolevs, beta_sobolevs):
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


def calculate_beta_sobolev(tau_sobolevs):
    """Calculate the beta Sobolev values based on the provided tau_sobolevs.
    Values from the previous iteration can be provided.

    Parameters
    ----------
    tau_sobolevs : pd.DataFrame
        Tau Sobolev opacities.

    Returns
    -------
    pd.DataFrame
        The latest Beta Sobolev opacities.
    """
    beta_sobolev = pd.DataFrame(
        0.0, index=tau_sobolevs.index, columns=tau_sobolevs.columns
    )

    numba_calculate_beta_sobolev(
        tau_sobolevs.values.ravel(), beta_sobolev.values.ravel()
    )

    return beta_sobolev


class TauSobolev(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    tau_sobolev : Pandas DataFrame, dtype float
          Sobolev optical depth for each line. Indexed by line.
          Columns as zones.
    """

    outputs = ("tau_sobolevs",)
    latex_name = (r"\tau_{\textrm{sobolev}}",)
    latex_formula = (
        r"\dfrac{\pi e^{2}}{m_{e} c}f_{lu}\lambda t_{exp}\
        n_{lower} \Big(1-\dfrac{g_{lower}n_{upper}}{g_{upper}n_{lower}}\Big)",
    )

    def calculate(
        self,
        lines,
        level_number_density,
        time_explosion,
        stimulated_emission_factor,
    ):
        """
        Calculate Sobolev line opacity.

        Calculates the Sobolev line opacity based on the provided parameters.

        Parameters
        ----------
        lines : pandas.DataFrame
            DataFrame containing information about spectral lines.
        level_number_density : pandas.DataFrame
            DataFrame with level number densities.
        time_explosion : astropy.units.Quantity
            Time since explosion.
        stimulated_emission_factor : float
            Factor for stimulated emission.

        Returns
        -------
        pandas.DataFrame
            Calculated Sobolev line opacity values.

        Raises
        ------
        ValueError
            If any calculated tau_sobolevs are nan or inf.
        """
        return calculate_sobolev_line_opacity(
            lines,
            level_number_density,
            time_explosion,
            stimulated_emission_factor,
        )


class BetaSobolev(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    beta_sobolev : Numpy Array, dtype float
    """

    outputs = ("beta_sobolev",)
    latex_name = (r"\beta_{\textrm{sobolev}}",)

    def calculate(self, tau_sobolevs):
        return calculate_beta_sobolev(tau_sobolevs)
