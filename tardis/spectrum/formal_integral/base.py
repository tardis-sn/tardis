import numpy as np
from scipy.interpolate import interp1d
import scipy.sparse as sp
import scipy.sparse.linalg as linalg

import warnings

from tardis import constants as const
from tardis.transport.montecarlo.configuration import montecarlo_globals

C_INV = 3.33564e-11
KB_CGS = 1.3806488e-16
H_CGS = 6.62606957e-27


# here will go function that are independent of the type of formal integral
class BoundsError(IndexError):
    pass


class IntegrationError(Exception):
    pass


def check(simulation_state, plasma, transport, raises=True):
    """
    A method that determines if the formal integral can be performed with
    the current configuration settings

    The function returns False if the configuration conflicts with the
    required settings. If raises evaluates to True, then a
    IntegrationError is raised instead
    """

    def raise_or_return(message):
        if raises:
            raise IntegrationError(message)
        else:
            warnings.warn(message)
            return False

    for obj in (simulation_state, plasma, transport):
        if obj is None:
            return raise_or_return(
                "The integrator is missing either model, plasma or "
                "transport. Please make sure these are provided to the "
                "FormalIntegrator."
            )

    if transport.line_interaction_type not in [
        "downbranch",
        "macroatom",
    ]:
        return raise_or_return(
            "The FormalIntegrator currently only works for "
            'line_interaction_type == "downbranch"'
            'and line_interaction_type == "macroatom"'
        )

    if montecarlo_globals.CONTINUUM_PROCESSES_ENABLED:
        return raise_or_return(
            "The FormalIntegrator currently does not work for "
            "continuum interactions."
        )

    return True


def calculate_p_values(R_max, N):
    """
    Calculates the p values of N

    Parameters
    ----------
    R_max : float64
    N : int64

    Returns
    -------
    float64
    """
    return np.arange(N).astype(np.float64) * R_max / (N - 1)


def intensity_black_body(nu, temperature):
    """
    Calculate the blackbody intensity.

    Parameters
    ----------
    nu : float64
        frequency
    temperature : float64
        Temperature

    Returns
    -------
    float64
    """
    if nu == 0:
        return np.nan  # to avoid ZeroDivisionError
    beta_rad = 1 / (KB_CGS * temperature)
    coefficient = 2 * H_CGS * C_INV * C_INV
    return coefficient * nu * nu * nu / (np.exp(H_CGS * nu * beta_rad) - 1)


def interpolate_integrator_quantities(
    source_function_state,
    interpolate_shells,
    simulation_state,
    transport,
    opacity_state,
    electron_densities,
):
    """Interpolate the integrator quantities to interpolate_shells.

    Parameters
    ----------
    source_function_state: SourceFunctionState
        Data class to hold the computed source function values
    interpolate_shells : int
        number of shells to interpolate to
    simulation_state : tardis.model.SimulationState
    transport : tardis.transport.montecarlo.MonteCarloTransportSolver
    opacity_state : OpacityStateNumba
    electron_densities : pd.DataFrame

    Returns
    -------
    tuple
        Interpolated values of att_S_ul, Jred_lu, Jbluelu, and e_dot_u
    """

    mct_state = transport.transport_state

    att_S_ul, Jred_lu, Jblue_lu, e_dot_u = (
            source_function_state.att_S_ul,
            source_function_state.Jred_lu,
            source_function_state.Jblue_lu,
            source_function_state.e_dot_u,
        )

    nshells = interpolate_shells
    r_middle = (
        mct_state.geometry_state.r_inner + mct_state.geometry_state.r_outer
    ) / 2.0

    r_integ = np.linspace(
        mct_state.geometry_state.r_inner[0],
        mct_state.geometry_state.r_outer[-1],
        nshells,
    )
    transport.r_inner_i = r_integ[:-1]
    transport.r_outer_i = r_integ[1:]

    r_middle_integ = (r_integ[:-1] + r_integ[1:]) / 2.0

    transport.electron_densities_integ = interp1d(
        r_middle,
        electron_densities.iloc[
            simulation_state.geometry.v_inner_boundary_index : simulation_state.geometry.v_outer_boundary_index
        ],
        fill_value="extrapolate",
        kind="nearest",
    )(r_middle_integ)
    # Assume tau_sobolevs to be constant within a shell
    # (as in the MC simulation)
    transport.tau_sobolevs_integ = interp1d(
        r_middle,
        opacity_state.tau_sobolev[
            :,
            simulation_state.geometry.v_inner_boundary_index : simulation_state.geometry.v_outer_boundary_index,
        ],
        fill_value="extrapolate",
        kind="nearest",
    )(r_middle_integ)
    att_S_ul = interp1d(r_middle, att_S_ul, fill_value="extrapolate")(
        r_middle_integ
    )
    Jred_lu = interp1d(r_middle, Jred_lu, fill_value="extrapolate")(
        r_middle_integ
    )
    Jblue_lu = interp1d(r_middle, Jblue_lu, fill_value="extrapolate")(
        r_middle_integ
    )
    e_dot_u = interp1d(r_middle, e_dot_u, fill_value="extrapolate")(
        r_middle_integ
    )

    # Set negative values from the extrapolation to zero
    att_S_ul = att_S_ul.clip(0.0)
    Jblue_lu = Jblue_lu.clip(0.0)
    Jred_lu = Jred_lu.clip(0.0)
    e_dot_u = e_dot_u.clip(0.0)
    return att_S_ul, Jred_lu, Jblue_lu, e_dot_u