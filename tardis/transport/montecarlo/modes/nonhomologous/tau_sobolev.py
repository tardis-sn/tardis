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
    velocity_gradient,
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
    velocity_gradient : np.ndarray
        Local dv/dr per cell.
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
    sobolev_optical_depth_coefficient = (
        (lines.wavelength_cm * lines.f_lu).values[np.newaxis].T
        * SOBOLEV_COEFFICIENT
        * stimulated_emission_factor
        * level_number_density.reindex(lines.droplevel(-1).index).values
    )
    velocity_gradient = np.abs(velocity_gradient.to(1 / u.s).value)
    tau_sobolevs = np.divide(
        sobolev_optical_depth_coefficient,
        velocity_gradient,
        out=np.full_like(sobolev_optical_depth_coefficient, np.inf),
        where=velocity_gradient != 0.0,
    )
    tau_sobolevs[sobolev_optical_depth_coefficient == 0.0] = 0.0

    if np.any(np.isnan(tau_sobolevs)):
        raise ValueError(
            "Some tau_sobolevs are nan in tau_sobolevs."
            " Something went wrong!"
        )

    return pd.DataFrame(
        tau_sobolevs,
        index=lines.index,
        columns=np.array(level_number_density.columns),
    )


def calculate_sobolev_optical_depth_coefficient(
    lines: pd.DataFrame,
    level_number_density: pd.DataFrame,
    stimulated_emission_factor: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate the velocity-gradient-independent Sobolev optical-depth coefficient.

    Parameters
    ----------
    lines : pandas.DataFrame
        Atomic line data containing wavelength and oscillator strength.
    level_number_density : pandas.DataFrame
        Lower-level number densities for each shell [cm^-3].
    stimulated_emission_factor : pandas.DataFrame
        Stimulated-emission correction for each line and shell.

    Returns
    -------
    pandas.DataFrame
        Coefficient for each line and shell [s^-1]. Dividing by the absolute
        projected velocity gradient gives the directional Sobolev optical
        depth.
    """
    sobolev_optical_depth_coefficient = (
        (lines.wavelength_cm * lines.f_lu).values[np.newaxis].T
        * SOBOLEV_COEFFICIENT
        * stimulated_emission_factor
        * level_number_density.reindex(lines.droplevel(-1).index).values
    )

    if np.any(np.isnan(sobolev_optical_depth_coefficient)) or np.any(
        np.isinf(np.abs(sobolev_optical_depth_coefficient))
    ):
        raise ValueError(
            "Some Sobolev optical-depth coefficients are nan, inf, or -inf."
        )

    return pd.DataFrame(
        sobolev_optical_depth_coefficient,
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
    tau_values_1d = tau_sobolevs.to_numpy().ravel()
    beta_values_1d = np.zeros(tau_values_1d.shape[0], dtype=np.float64)

    numba_calculate_beta_sobolev(tau_values_1d, beta_values_1d)

    return pd.DataFrame(
        beta_values_1d.reshape(tau_sobolevs.shape),
        index=tau_sobolevs.index,
        columns=tau_sobolevs.columns,
    )


def calculate_beta_sobolev_directional(
    sobolev_optical_depth_coefficient: pd.DataFrame,
    velocity_gradient: u.Quantity,
    velocity_over_radius: u.Quantity,
) -> pd.DataFrame:
    """Calculate the angle-averaged Sobolev escape probability.

    Parameters
    ----------
    sobolev_optical_depth_coefficient : pandas.DataFrame
        Velocity-gradient-independent Sobolev optical-depth coefficient [s^-1].
    velocity_gradient : astropy.units.Quantity
        Radial velocity derivative in each shell [s^-1].
    velocity_over_radius : astropy.units.Quantity
        Ratio of radial velocity to radius in each shell [s^-1].

    Returns
    -------
    pandas.DataFrame
        Angle-averaged Sobolev escape probabilities for each line and shell.
    """
    mus, weights = np.polynomial.legendre.leggauss(20)
    beta_sobolevs = np.zeros_like(
        sobolev_optical_depth_coefficient.values
    )
    velocity_gradient = velocity_gradient.to(1 / u.s).value
    velocity_over_radius = velocity_over_radius.to(1 / u.s).value

    for mu, weight in zip(mus, weights, strict=True):
        projected_velocity_gradient = (
            mu * mu * velocity_gradient
            + (1.0 - mu * mu) * velocity_over_radius
        )
        tau_sobolevs = np.divide(
            sobolev_optical_depth_coefficient.values,
            np.abs(projected_velocity_gradient),
            out=np.full_like(
                sobolev_optical_depth_coefficient.values, np.inf
            ),
            where=projected_velocity_gradient != 0.0,
        )
        tau_sobolevs[
            sobolev_optical_depth_coefficient.values == 0.0
        ] = 0.0

        beta_directional = np.empty_like(tau_sobolevs)
        thick = tau_sobolevs > 1e3
        thin = tau_sobolevs < 1e-4
        intermediate = np.logical_not(np.logical_or(thick, thin))
        beta_directional[thick] = tau_sobolevs[thick] ** -1
        beta_directional[thin] = 1.0 - 0.5 * tau_sobolevs[thin]
        beta_directional[intermediate] = (
            1.0 - np.exp(-tau_sobolevs[intermediate])
        ) / tau_sobolevs[intermediate]
        beta_sobolevs += 0.5 * weight * beta_directional

    return pd.DataFrame(
        beta_sobolevs,
        index=sobolev_optical_depth_coefficient.index,
        columns=sobolev_optical_depth_coefficient.columns,
    )


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
        velocity_gradient,
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
        velocity_gradient : np.ndarray
            Local dv/dr per cell.
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
            velocity_gradient,
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
