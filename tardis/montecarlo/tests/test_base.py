import os
import pandas as pd
import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_almost_equal

###
# Save and Load
###


@pytest.fixture(scope="module", autouse=True)
def to_hdf_buffer(hdf_file_path, simulation_verysimple):
    simulation_verysimple.transport.to_hdf(
        hdf_file_path, name="transport", overwrite=True
    )


transport_properties = [
    "output_nu",
    "output_energy",
    "nu_bar_estimator",
    "j_estimator",
    "montecarlo_virtual_luminosity",
    "last_interaction_in_nu",
    "last_interaction_type",
    "last_line_interaction_in_id",
    "last_line_interaction_in_shell_id",
    "last_line_interaction_out_id",
    "last_line_interaction_out_shell_id",
    "packet_luminosity",
]


@pytest.mark.parametrize("attr", transport_properties)
def test_hdf_transport(hdf_file_path, simulation_verysimple, attr):
    actual = getattr(simulation_verysimple.transport, attr)
    if hasattr(actual, "cgs"):
        actual = actual.cgs.value
    path = os.path.join("transport", attr)
    expected = pd.read_hdf(hdf_file_path, path)
    assert_almost_equal(actual, expected.values)
