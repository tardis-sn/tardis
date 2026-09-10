# Utility functions for continuum calculations

import pandas as pd

__all__ = ["get_ion_multi_index"]


def get_ion_multi_index(multi_index_full, next_higher=True):
    atomic_number = multi_index_full.get_level_values(0)
    ion_number = multi_index_full.get_level_values(1)
    if next_higher is True:
        ion_number += 1
    return pd.MultiIndex.from_arrays([atomic_number, ion_number])
