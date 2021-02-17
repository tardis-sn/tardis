
import os
import yaml
import pytest

from tardis.io.atom_data.base import AtomData
from tardis.simulation import Simulation
from tardis.io.config_reader import Configuration

from tardis.montecarlo.montecarlo_numba.base import montecarlo_main_loop as main_loop
from tardis.montecarlo.montecarlo_numba.base import montecarlo_radial1d
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    # PacketCollection,
    # VPacketCollection,
    NumbaModel,
    numba_plasma_initialize,
    # Estimators,
    # configuration_initialize,
)

# @pytest.mark.skipif(
#     'not config.getvalue("integration-tests")',
#     reason="integration tests are not included in this run",
# )
# @pytest.mark.integration
# class TestIntegration(object):
#     """Slow integration test for various setups present in subdirectories of
#     ``tardis/tests/integration_tests``.
#     """

#@classmethod
# @pytest.fixture(scope="class")
def test_setup(data_path):
    """
    This method does initial setup of creating configuration and performing
    a single run of integration test.
    """
    # Get the path to HDF file:
    config_file = os.path.join(
        data_path["config_dirpath"], "config.yml"
    )
    # Load atom data file separately, pass it for forming tardis config.

    atom_data_name = yaml.load(open(config_file), Loader=yaml.CLoader)[
        "atom_data"
    ]
    atom_data_filepath = '/Users/ghost/Desktop/tardis/tardis/reference-data-full/tardis-refdata/atom_data/kurucz_cd23_chianti_H_He.h5'
    print(os.path.exists('/Users/ghost/Desktop/tardis/tardis/reference-data-full/tardis-refdata/atom_data/kurucz_cd23_chianti_H_He.h5'), atom_data_filepath)
    atom_data = AtomData.from_hdf(atom_data_filepath)
    tardis_config = Configuration.from_yaml(config_file)

    sim = Simulation.from_config(tardis_config, atom_data=atom_data)

    # Init model
    numba_plasma = numba_plasma_initialize(
        sim.plasma, line_interaction_type="macroatom"
    )

    runner = sim.runner
    model = sim.model
    # print("electron density values: ", model.electron_densities.values)
    # pass numba_plasma, runner to montecarlo_radial1d
    runner._initialize_geometry_arrays(model)
    runner._initialize_estimator_arrays(numba_plasma.tau_sobolev.shape)
    runner._initialize_packets(model.t_inner.value, 100000, 0)
    montecarlo_radial1d(model, numba_plasma, runner)

    print("runner values", dir(runner))
    assert 1 == 1