import pytest
import pandas as pd
import os
import numpy.testing as npt
from copy import deepcopy
from tardis.base import run_tardis
from pandas.testing import assert_frame_equal

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
from tardis.simulation import Simulation


@pytest.mark.xfail(reason="To be implemented")
def test_montecarlo_radial1d():
    assert False


def test_montecarlo_main_loop(
    config_montecarlo_1e5_verysimple,
    atomic_dataset,
    tardis_ref_path,
    tmpdir,
    set_seed_fixture,
    random_call_fixture,
    request,
):

    montecarlo_configuration.LEGACY_MODE_ENABLED = True
    # Setup model config from verysimple
    atomic_data = deepcopy(atomic_dataset)
    config_montecarlo_1e5_verysimple.montecarlo.last_no_of_packets = 1e5
    config_montecarlo_1e5_verysimple.montecarlo.no_of_virtual_packets = 0
    config_montecarlo_1e5_verysimple.montecarlo.iterations = 1
    config_montecarlo_1e5_verysimple.plasma.line_interaction_type = "macroatom"
    del config_montecarlo_1e5_verysimple["config_dirname"]

    sim = Simulation.from_config(
        config_montecarlo_1e5_verysimple, atom_data=atomic_data
    )
    sim.run()

    compare_fname = os.path.join(
        tardis_ref_path, "montecarlo_1e5_compare_data.h5"
    )
    if request.config.getoption("--generate-reference"):
        sim.to_hdf(compare_fname, overwrite=True)

    # Load compare data from refdata
    expected_nu = pd.read_hdf(
        compare_fname, key="/simulation/runner/output_nu"
    ).values
    expected_energy = pd.read_hdf(
        compare_fname, key="/simulation/runner/output_energy"
    ).values
    expected_nu_bar_estimator = pd.read_hdf(
        compare_fname, key="/simulation/runner/nu_bar_estimator"
    ).values
    expected_j_estimator = pd.read_hdf(
        compare_fname, key="/simulation/runner/j_estimator"
    ).values

    actual_energy = sim.runner.output_energy
    actual_nu = sim.runner.output_nu
    actual_nu_bar_estimator = sim.runner.nu_bar_estimator
    actual_j_estimator = sim.runner.j_estimator

    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-13
    )
    npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-13)
    npt.assert_allclose(actual_energy.value, expected_energy, rtol=1e-13)
    npt.assert_allclose(actual_nu.value, expected_nu, rtol=1e-13)
