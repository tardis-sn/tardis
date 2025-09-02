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


def check(simulation_state, opacity_state, transport, raises=True):
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

    for obj in (simulation_state, opacity_state, transport):
        if obj is None:
            return raise_or_return(
                "The integrator is missing either model, opacity state or "
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
