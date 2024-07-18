import os
import pytest
from tardis.io.configuration.config_reader import Configuration
from tardis.simulation import Simulation
import astropy.units as u

if "QT_API" in os.environ:
    from PyQt5 import QtWidgets
    from tardis.gui.widgets import Tardis
    from tardis.gui.datahandler import SimpleTableModel


@pytest.fixture(scope="module")
def config():
    return Configuration.from_yaml(
        "tardis/io/tests/data/tardis_configv1_verysimple.yml"
    )


@pytest.fixture(scope="module")
def simulation_one_loop(
    atomic_data_fname, config, tardis_ref_data, generate_reference
):
    config.atom_data = atomic_data_fname
    config.montecarlo.iterations = 2
    config.montecarlo.no_of_packets = int(4e4)
    config.montecarlo.last_no_of_packets = int(4e4)

    simulation = Simulation.from_config(config)
    simulation.run_convergence()
    simulation.run_final()

    return simulation


@pytest.mark.skipif(
    "QT_API" not in os.environ, reason="enviroment variable QT_API is not set"
)
def test_gui(simulation_one_loop):
    simulation = simulation_one_loop
    app = QtWidgets.QApplication([])
    tablemodel = SimpleTableModel
    win = Tardis(tablemodel)
    win.show_model(simulation)
    app.quit()
