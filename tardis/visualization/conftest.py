from tardis.base import run_tardis
import pytest
from copy import deepcopy
from tardis.base import run_tardis
from tardis.workflows.standard_tardis_workflow import StandardTARDISWorkflow

@pytest.fixture(scope="module")
def simulation_simple_tracked(config_verysimple, atomic_dataset):
    """
    Instantiate SDEC plotter using a simple simulation model.

    Parameters
    ----------
    config_verysimple : tardis.io.config_reader.Configuration
        Configuration object for a very simple simulation.
    atomic_dataset : str or tardis.atomic.AtomData
        Atomic data.

    Returns
    -------
    sim: tardis.simulation.base.Simulation
        Simulation object.
    """
    # Setup simulation configuration using config_verysimple and
    # override properties in such a way to make the simulation run faster
    config_verysimple.montecarlo.iterations = 3
    config_verysimple.montecarlo.no_of_packets = 4000
    config_verysimple.montecarlo.last_no_of_packets = -1
    config_verysimple.spectrum.virtual.virtual_packet_logging = True
    config_verysimple.montecarlo.no_of_virtual_packets = 1
    config_verysimple.spectrum.num = 2000
    config_verysimple.montecarlo.tracking.track_rpacket = True
    atomic_data = deepcopy(atomic_dataset)
    sim = run_tardis(
        config_verysimple,
        atom_data=atomic_data,
        show_convergence_plots=False,
        log_level="CRITICAl",
    )
    return sim

@pytest.fixture(scope="class")
def workflow_simple_tracked(config_verysimple, atomic_data_fname):
    config = deepcopy(config_verysimple)
    config.atom_data = atomic_data_fname
    config.montecarlo.iterations = 3
    config.montecarlo.no_of_packets = 4000
    config.montecarlo.last_no_of_packets = -1
    config.spectrum.virtual.virtual_packet_logging = True
    config.montecarlo.no_of_virtual_packets = 1
    config.spectrum.num = 2000
    config_verysimple.montecarlo.tracking.track_rpacket = True
    
    workflow = StandardTARDISWorkflow(config, enable_virtual_packet_logging=True)
    workflow.run()
    return workflow
