import numpy as np
import pandas as pd
import os

def data_path():
    return os.path.join('tardis', 'plasma', 'tests', 'data', 'ref_data.hdf')


def plasma_compare_data(path):
    fpath = data_path()
    with pd.HDFStore(fpath, 'r') as plasma:
        data = plasma[path]
    return data

def test_level_boltzmann_factor_lte(level_boltzmann_factor_lte, levels):
    expected = plasma_compare_data('/plasma/level_boltzmann_factor_lte')
    assert np.allclose(level_boltzmann_factor_lte.ix[2].ix[0].ix[0], expected.ix[2].ix[0].ix[0])
    assert np.allclose(level_boltzmann_factor_lte.ix[2].ix[1].ix[0], expected.ix[2].ix[1].ix[0])
    assert np.allclose(level_boltzmann_factor_lte.ix[2].ix[0].ix[10],
        expected.ix[2].ix[0].ix[10])

def test_level_boltzmann_factor_dilute_lte(level_boltzmann_factor_dilute_lte,
    levels):
    expected = plasma_compare_data('/plasma/level_boltzmann_factor_dilute_lte')
    assert np.allclose(level_boltzmann_factor_dilute_lte.ix[2].ix[0].ix[0], expected.ix[2].ix[0].ix[0])
    assert np.allclose(level_boltzmann_factor_dilute_lte.ix[2].ix[1].ix[0], expected.ix[2].ix[1].ix[0])
    assert np.allclose(level_boltzmann_factor_dilute_lte.ix[2].ix[0].ix[10],
        expected.ix[2].ix[0].ix[10])

def test_lte_partition_function(partition_function, levels):
    expected = plasma_compare_data('/plasma/partition_function')
    assert np.allclose(partition_function.ix[2].ix[0], expected.ix[2].ix[0])
    assert np.allclose(partition_function.ix[2].ix[1], expected.ix[2].ix[1])
    assert np.allclose(partition_function.ix[2].ix[2], expected.ix[2].ix[2])

