from copy import deepcopy

import numpy as np
import pytest
from numba import njit

from tardis.opacities.opacity_state import opacity_state_initialize
from tardis.simulation import Simulation
from tardis.transport.montecarlo import RPacket
from tardis.transport.montecarlo.estimators.radfield_mc_estimators import (
    RadiationFieldMCEstimators,
)
from tardis.opacities.opacity_state_numba import (
    opacity_state_numba_initialize,
)
from tardis.transport.montecarlo.packet_collections import (
    VPacketCollection,
)
from tardis.transport.montecarlo.weighted_packet_source import (
    BlackBodyWeightedSource,
)


@pytest.fixture(scope="function")
def montecarlo_main_loop_config(
    config_montecarlo_1e5_verysimple,
):
    # Setup model config from verysimple

    config_montecarlo_1e5_verysimple.montecarlo.last_no_of_packets = 1e5
    config_montecarlo_1e5_verysimple.montecarlo.no_of_virtual_packets = 0
    config_montecarlo_1e5_verysimple.montecarlo.iterations = 1
    config_montecarlo_1e5_verysimple.plasma.line_interaction_type = "macroatom"

    del config_montecarlo_1e5_verysimple["config_dirname"]
    return config_montecarlo_1e5_verysimple


@pytest.fixture(scope="package")
def nb_simulation_verysimple(config_verysimple, atomic_dataset):
    atomic_data = deepcopy(atomic_dataset)
    sim = Simulation.from_config(config_verysimple, atom_data=atomic_data)
    sim.iterate(10)
    return sim


@pytest.fixture(scope="package")
def simple_weighted_packet_source():
    packet_source = BlackBodyWeightedSource(base_seed=1963)

    return packet_source


@pytest.fixture(scope="package")
def verysimple_opacity_state(nb_simulation_verysimple):
    return opacity_state_numba_initialize(
        nb_simulation_verysimple.plasma,
        line_interaction_type="macroatom",
        disable_line_scattering=False,
    )


@pytest.fixture(scope="package")
def verysimple_numba_radial_1d_geometry(nb_simulation_verysimple):
    return nb_simulation_verysimple.simulation_state.geometry.to_numba()


@pytest.fixture(scope="package")
def verysimple_time_explosion(nb_simulation_verysimple):
    model = nb_simulation_verysimple.simulation_state
    return model.time_explosion.cgs.value


@pytest.fixture(scope="package")
def verysimple_estimators(nb_simulation_verysimple):
    transport = nb_simulation_verysimple.transport

    return RadiationFieldMCEstimators(
        transport.j_estimator,
        transport.nu_bar_estimator,
        transport.j_blue_estimator,
        transport.Edotlu_estimator,
        transport.photo_ion_estimator,
        transport.stim_recomb_estimator,
        transport.bf_heating_estimator,
        transport.stim_recomb_cooling_estimator,
        transport.photo_ion_estimator_statistics,
    )


@pytest.fixture(scope="package")
def verysimple_vpacket_collection(nb_simulation_verysimple):
    spectrum_frequency_grid = (
        nb_simulation_verysimple.transport.spectrum_frequency_grid.value
    )
    return VPacketCollection(
        source_rpacket_index=0,
        spectrum_frequency_grid=spectrum_frequency_grid,
        number_of_vpackets=0,
        v_packet_spawn_start_frequency=0,
        v_packet_spawn_end_frequency=np.inf,
        temporary_v_packet_bins=0,
    )


@pytest.fixture(scope="package")
def verysimple_3vpacket_collection(nb_simulation_verysimple):
    spectrum_frequency_grid = (
        nb_simulation_verysimple.transport.spectrum_frequency_grid.value
    )
    return VPacketCollection(
        source_rpacket_index=0,
        spectrum_frequency_grid=spectrum_frequency_grid,
        number_of_vpackets=3,
        v_packet_spawn_start_frequency=0,
        v_packet_spawn_end_frequency=np.inf,
        temporary_v_packet_bins=0,
    )


@pytest.fixture(scope="package")
def verysimple_packet_collection(nb_simulation_verysimple):
    return nb_simulation_verysimple.transport.transport_state.packet_collection


@pytest.fixture(scope="function")
def packet(verysimple_packet_collection):
    return RPacket(
        r=7.5e14,
        nu=verysimple_packet_collection.initial_nus[0],
        mu=verysimple_packet_collection.initial_mus[0],
        energy=verysimple_packet_collection.initial_energies[0],
        seed=1963,
        index=0,
    )


@pytest.fixture(scope="function")
def static_packet():
    return RPacket(
        r=7.5e14,
        nu=0.4,
        mu=0.3,
        energy=0.9,
        seed=1963,
        index=0,
    )


@pytest.fixture
def set_seed_fixture():
    def set_seed(value):
        np.random.seed(value)

    return njit(set_seed)


@pytest.fixture
def random_call_fixture():
    def random_call():
        np.random.random()

    return njit(random_call)
