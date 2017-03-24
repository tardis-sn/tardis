import os
import pytest
import h5py
import numpy as np

from tardis import __path__ as path

@pytest.fixture(scope='module')
def partition_compare_data_fname():
    fname = 'partition_compare_data.h5'
    return os.path.join(path[0], 'plasma', 'tests', 'data', fname)

@pytest.fixture()
def partition_compare_data(partition_compare_data_fname):
    return h5py.File(partition_compare_data_fname, 'r')

def test_level_boltzmann_factor_lte(level_boltzmann_factor_lte, partition_compare_data, levels):
    expected_lbfl = partition_compare_data['lbfl']
    assert np.allclose(level_boltzmann_factor_lte.ix[2].ix[0].ix[0], expected_lbfl[0])
    assert np.allclose(level_boltzmann_factor_lte.ix[2].ix[1].ix[0], expected_lbfl[1])
    assert np.allclose(level_boltzmann_factor_lte.ix[2].ix[0].ix[10], expected_lbfl[2])

def test_level_boltzmann_factor_dilute_lte(level_boltzmann_factor_dilute_lte, partition_compare_data,
    levels):
    expected_lbfdl = partition_compare_data['lbfdl']
    assert np.allclose(level_boltzmann_factor_dilute_lte.ix[2].ix[0].ix[0], expected_lbfdl[0])
    assert np.allclose(level_boltzmann_factor_dilute_lte.ix[2].ix[1].ix[0], expected_lbfdl[1])
    assert np.allclose(level_boltzmann_factor_dilute_lte.ix[2].ix[0].ix[10], expected_lbfdl[2])

def test_lte_partition_function(partition_function, partition_compare_data, levels):
    expected_pf = partition_compare_data['pf']
    assert np.allclose(partition_function.ix[2].ix[0], expected_pf[0])
    assert np.allclose(partition_function.ix[2].ix[1], expected_pf[1])
    assert np.allclose(partition_function.ix[2].ix[2], expected_pf[2])



