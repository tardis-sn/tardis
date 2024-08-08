from dataclasses import dataclass

import pandas as pd


@dataclass
class ChiantiCollisionData:
    """Collisional data from Chianti sources
    data : (pandas.DataFrame, np.array)
    A DataFrame containing the *electron collisions data* with:
        index : atomic_number, ion_number, level_number_lower, level_number_upper
        columns : e_col_id, delta_e, g_ratio, c_ul;
    temperatures : np.array
        An array with the collision temperatures.
    """

    data: pd.DataFrame
    temperatures: pd.DataFrame


@dataclass
class CMFGENCollisionData:
    """Collisional data from CMFGEN sources
    data : (pandas.DataFrame, np.array)
        yg data from CMFGEN
    temperatures : np.array
        An array with the collision temperatures.
    """

    data: pd.DataFrame
    temperatures: pd.DataFrame
