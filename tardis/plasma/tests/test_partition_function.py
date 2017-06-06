import numpy as np
import pandas as pd
import os
import pytest

ref_data_path=""

@pytest.fixture(autouse=True)
def set_data_path(compare_with_reference):
    global ref_data_path
    if compare_with_reference:
        ref_data_path= compare_with_reference
    else:
        ref_data_path= os.path.join('tardis', 'plasma', 'tests', 'data', 'ref_data.h5')


def plasma_compare_data(path):
    with pd.HDFStore(ref_data_path, 'r') as plasma:
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

