import os
import pandas as pd
import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_almost_equal
from pathlib import Path

###
# Save and Load
###


@pytest.fixture(scope="module", autouse=True)
def to_hdf_buffer(hdf_file_path, simulation_verysimple):
    simulation_verysimple.transport.to_hdf(
        hdf_file_path, name="transport", overwrite=True
    )
    simulation_verysimple.transport.transport_state.to_hdf(
        hdf_file_path, name="transport_state", overwrite=True
    )


transport_properties = [
    "last_interaction_in_nu",
    "last_interaction_type",
    "last_line_interaction_in_id",
    "last_line_interaction_out_id",
    "last_line_interaction_shell_id",
]


@pytest.mark.parametrize("attr", transport_properties)
def test_hdf_transport(hdf_file_path, simulation_verysimple, attr):
    actual = getattr(simulation_verysimple.transport, attr)
    if hasattr(actual, "cgs"):
        actual = actual.cgs.value
    path = f"transport/{attr}"
    expected = pd.read_hdf(hdf_file_path, path)
    assert_almost_equal(actual, expected.values)


transport_state_properties = [
    "output_nu",
    "output_energy",
    "nu_bar_estimator",
    "j_estimator",
    "montecarlo_virtual_luminosity",
    "packet_luminosity",
]


@pytest.mark.parametrize("attr", transport_state_properties)
def test_hdf_transport_state(hdf_file_path, simulation_verysimple, attr):
    actual = getattr(simulation_verysimple.transport.transport_state, attr)
    if hasattr(actual, "cgs"):
        actual = actual.cgs.value
    path = f"transport_state/{attr}"
    expected = pd.read_hdf(hdf_file_path, path)
    assert_almost_equal(actual, expected.values)
