from copy import deepcopy

import pytest

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
    # Use deepcopy to avoid mutating the session-scoped config_verysimple
    config = deepcopy(config_verysimple)
    config.montecarlo.iterations = 3
    config.montecarlo.no_of_packets = 4000
    config.montecarlo.last_no_of_packets = -1
    # Default to vpackets=0 for performance
    config.spectrum.virtual.virtual_packet_logging = False
    config.montecarlo.no_of_virtual_packets = 0
    config.spectrum.num = 2000
    config.montecarlo.tracking.track_rpacket = True
    atomic_data = deepcopy(atomic_dataset)
    sim = run_tardis(
        config,
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
    # Default to vpackets=0 for performance
    config.spectrum.virtual.virtual_packet_logging = False
    config.montecarlo.no_of_virtual_packets = 0
    config.spectrum.num = 2000
    config.montecarlo.tracking.track_rpacket = True

    workflow = StandardTARDISWorkflow(
        config, enable_virtual_packet_logging=False
    )
    workflow.run()
    return workflow


@pytest.fixture(scope="module")
def simulation_simple_tracked_with_vpackets(config_verysimple, atomic_dataset):
    """
    Instantiate SDEC plotter using a simple simulation model with vpackets enabled.

    Parameters
    ----------
    config_verysimple : tardis.io.config_reader.Configuration
        Configuration object for a very simple simulation.
    atomic_dataset : str or tardis.atomic.AtomData
        Atomic data.

    Returns
    -------
    sim: tardis.simulation.base.Simulation
        Simulation object with virtual packets enabled.
    """
    # Use deepcopy to avoid mutating the session-scoped config_verysimple
    config = deepcopy(config_verysimple)
    config.montecarlo.iterations = 3
    config.montecarlo.no_of_packets = 4000
    config.montecarlo.last_no_of_packets = -1
    # Explicitly enable vpackets for vpacket-specific tests
    config.spectrum.virtual.virtual_packet_logging = True
    config.montecarlo.no_of_virtual_packets = 5
    config.spectrum.num = 2000
    config.montecarlo.tracking.track_rpacket = True
    atomic_data = deepcopy(atomic_dataset)
    sim = run_tardis(
        config,
        atom_data=atomic_data,
        show_convergence_plots=False,
        log_level="CRITICAl",
    )
    return sim


@pytest.fixture(scope="class")
def workflow_simple_tracked_with_vpackets(config_verysimple, atomic_data_fname):
    """Workflow fixture with virtual packets enabled for vpacket-specific tests."""
    config = deepcopy(config_verysimple)
    config.atom_data = atomic_data_fname
    config.montecarlo.iterations = 3
    config.montecarlo.no_of_packets = 4000
    config.montecarlo.last_no_of_packets = -1
    # Explicitly enable vpackets for vpacket-specific tests
    config.spectrum.virtual.virtual_packet_logging = True
    config.montecarlo.no_of_virtual_packets = 5
    config.spectrum.num = 2000
    config.montecarlo.tracking.track_rpacket = True

    workflow = StandardTARDISWorkflow(
        config, enable_virtual_packet_logging=True
    )
    workflow.run()
    return workflow
