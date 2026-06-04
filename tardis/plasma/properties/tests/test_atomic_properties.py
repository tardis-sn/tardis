import os

import numpy as np
import numpy.testing as npt

import pandas as pd
import pandas.testing as pdt

import pytest
from astropy import units as u

from tardis.plasma.properties.atomic import * 
from tardis.plasma.base import BasePlasma


photo_ionization_data = PhotoIonizationData(plasma_parent=BasePlasma)

def test_photoionization_data_calculate(atomic_dataset,regression_data):
    actual_photoionization_data = photo_ionization_data.calculate(atomic_dataset,[])[0]
    expected_photoionization_data = regression_data.sync_dataframe(actual_photoionization_data,key="photo_ionization_data")
    pdt.assert_frame_equal(actual_photoionization_data,expected_photoionization_data)

    actual_block_references = photo_ionization_data.calculate(atomic_dataset,[])[1]
    expected_block_references = regression_data.sync_ndarray(actual_block_references)
    npt.assert_array_equal(actual_block_references,expected_block_references)


    actual_photo_ion_index = photo_ionization_data.calculate(atomic_dataset,[])[2]
    actual_photo_ion_index_np = actual_photo_ion_index.to_numpy()
    expected_photo_ion_index = regression_data.sync_ndarray(actual_photo_ion_index_np)
    npt.assert_array_equal(actual_photo_ion_index_np,expected_photo_ion_index)

    actual_nu_i = photo_ionization_data.calculate(atomic_dataset,[])[3]
    actual_nu_i_np = actual_nu_i.to_numpy()
    expected_nu_i = regression_data.sync_ndarray(actual_nu_i_np)
    npt.assert_array_equal(actual_nu_i_np,expected_nu_i)

    actual_energy_i = photo_ionization_data.calculate(atomic_dataset,[])[4]
    actual_energy_i_np = actual_energy_i.to_numpy()
    expected_energy_i = regression_data.sync_ndarray(actual_energy_i_np)
    npt.assert_array_equal(actual_energy_i_np,expected_energy_i)

    actual_photo_ion_idx = photo_ionization_data.calculate(atomic_dataset,[])[5]
    expected_photo_ion_idx = regression_data.sync_dataframe(actual_photo_ion_idx,key="photo_ion_idx")
    pdt.assert_frame_equal(actual_photo_ion_idx,expected_photo_ion_idx)

    actual_level2continuum_edge_idx = photo_ionization_data.calculate(atomic_dataset,[])[6]
    actual_level2continuum_edge_idx_np = actual_level2continuum_edge_idx.to_numpy()
    expected_level2continuum_edge_idx = regression_data.sync_ndarray(actual_level2continuum_edge_idx_np)
    npt.assert_array_equal(actual_level2continuum_edge_idx_np,expected_level2continuum_edge_idx)

    actual_level_idxs2continuum_idx = photo_ionization_data.calculate(atomic_dataset,[])[7]
    expected_level_idxs2continuum_idx = regression_data.sync_dataframe(actual_level_idxs2continuum_idx,key="level_idxs2continuum_idx")
    pdt.assert_frame_equal(actual_level_idxs2continuum_idx,expected_level_idxs2continuum_idx)

    
