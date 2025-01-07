import logging
import numpy as np
from astropy import units as u

logger = logging.getLogger(__name__)

def validate_radiative_temperature(t_radiative):
    """
    Validates the radiative temperature to ensure it is physically reasonable.

    Parameters
    ----------
    t_radiative : Quantity
        The radiative temperature array.

    Raises
    ------
    ValueError
        If any radiative temperature is below 1000 K.
    """
    if np.any(t_radiative < 1000 * u.K):
        min_t_rad = t_radiative[np.argmin(t_radiative)]
        min_shell = np.argmin(t_radiative)
        logging.critical(
            "Radiative temperature is too low in some of the shells, temperatures below 1000K "
            f"(e.g., T_rad = {min_t_rad} in shell {min_shell} in your model) "
            "are not accurately handled by TARDIS."
        )
        raise ValueError(
            f"Radiative temperature below 1000 K detected: T_rad = {min_t_rad} in shell {min_shell}."
        )