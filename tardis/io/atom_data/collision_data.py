from dataclasses import dataclass

import pandas as pd


@dataclass
class ChiantiCollisionData:
    """
    Collisional data from Chianti sources.

    Parameters
    ----------
    data : pandas.DataFrame
        A DataFrame containing the *electron collisions data* with:

        - index: atomic_number, ion_number, level_number_lower, level_number_upper
        - columns: e_col_id, delta_e, g_ratio, c_ul

    temperatures : numpy.ndarray
        An array with the collision temperatures.
    """

    data: pd.DataFrame
    temperatures: pd.DataFrame


@dataclass
class CMFGENCollisionData:
    """
    Collisional data from CMFGEN sources.

    Parameters
    ----------
    data : pandas.DataFrame
        yg data from CMFGEN

    temperatures : numpy.ndarray
        An array with the collision temperatures.
    """

    data: pd.DataFrame
    temperatures: pd.DataFrame
