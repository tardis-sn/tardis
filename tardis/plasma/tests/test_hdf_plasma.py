import os
import pandas as pd
import pytest

from numpy.testing import assert_almost_equal
import pandas.testing as pdt
from tardis.plasma.properties import property_collections

###
# saving and loading of plasma properties in the HDF file
###


# @pytest.fixture(scope="module", autouse=True)
# def to_hdf_buffer(hdf_file_path, simulation_verysimple):
#     simulation_verysimple.plasma.to_hdf(hdf_file_path, overwrite=True)


plasma_properties_list = [
    "number_density",
    "beta_rad",
    "general_level_boltzmann_factor",
    "level_boltzmann_factor",
    "stimulated_emission_factor",
    "t_electrons",
    "wavelength_cm",
    "lines_lower_level_index",
    "ionization_data",
    "density",
    "atomic_mass",
    "level_number_density",
    "lines_upper_level_index",
    "nu",
    "beta_sobolev",
    "transition_probabilities",
    "phi",
    "electron_densities",
    "t_rad",
    "selected_atoms",
    "ion_number_density",
    "partition_function",
    "abundance",
    "g_electron",
    "g",
    "lines",
    "f_lu",
    "tau_sobolevs",
    "j_blues",
    "metastability",
    "w",
    "excitation_energy",
]


@pytest.mark.parametrize("attr", plasma_properties_list)
def test_hdf_plasma(simulation_verysimple, attr, snapshot_np):
    if hasattr(simulation_verysimple.plasma, attr):
        actual = getattr(simulation_verysimple.plasma, attr)
        if hasattr(actual, "cgs"):
            actual = actual.cgs.value
        assert snapshot_np == actual


def test_hdf_levels(simulation_verysimple, snapshot_pd):
    actual = getattr(simulation_verysimple.plasma, "levels")
    if hasattr(actual, "cgs"):
        actual = actual.cgs.value
    assert snapshot_pd == pd.DataFrame(actual)


scalars_list = ["time_explosion", "link_t_rad_t_electron"]


@pytest.mark.parametrize("attr", scalars_list)
def test_hdf_scalars(simulation_verysimple, attr, snapshot_np):
    actual = getattr(simulation_verysimple.plasma, attr)
    if hasattr(actual, "cgs"):
        actual = actual.cgs.value
    assert snapshot_np == actual


def test_hdf_helium_treatment(hdf_file_path, simulation_verysimple, snapshot):
    actual = getattr(simulation_verysimple.plasma, "helium_treatment")
    assert snapshot == actual



def test_atomic_data_uuid(hdf_file_path, simulation_verysimple, snapshot):
    actual = getattr(simulation_verysimple.plasma.atomic_data, "uuid1")
    assert snapshot == actual


# @pytest.fixture(scope="module", autouse=True)
# def to_hdf_collection_buffer(hdf_file_path, simulation_verysimple):
#     simulation_verysimple.plasma.to_hdf(
#         hdf_file_path,
#         name="collection",
#         collection=property_collections.basic_inputs,
#         overwrite=True,
#     )


collection_properties = ["t_rad", "w", "density"]


@pytest.mark.parametrize("attr", collection_properties)
def test_collection(hdf_file_path, simulation_verysimple, attr, snapshot_np):
    actual = getattr(simulation_verysimple.plasma, attr)
    if hasattr(actual, "cgs"):
        actual = actual.cgs.value
    assert snapshot_np == actual

