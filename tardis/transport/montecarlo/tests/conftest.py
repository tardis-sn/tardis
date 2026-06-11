from copy import deepcopy

import numpy as np
import pytest

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.model.geometry.radial1d_nonhomologous import (
    NumbaNonhomologousRadial1DGeometry,
)
from tardis.opacities.opacity_state_numba import (
    OpacityStateNumba,
    opacity_state_numba_initialize,
)
from tardis.opacities.opacity_state_numba_iip import OpacityStateNumbaIIP
from tardis.simulation import Simulation
from tardis.transport.montecarlo import RPacket
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    init_estimators_bulk,
)
from tardis.transport.montecarlo.estimators.estimators_continuum import (
    init_estimators_continuum,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    VPacketCollection,
)


@pytest.fixture(scope="function")
def montecarlo_transport_config(
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
def verysimple_numba_nonhomologous_geometry(nb_simulation_verysimple):
    geometry = nb_simulation_verysimple.simulation_state.geometry
    return NumbaNonhomologousRadial1DGeometry(
        geometry.r_inner_active.to("cm").value,
        geometry.r_outer_active.to("cm").value,
        geometry.v_inner_active.to("cm/s").value,
        geometry.v_outer_active.to("cm/s").value,
    )


@pytest.fixture(scope="package")
def verysimple_time_explosion(nb_simulation_verysimple):
    model = nb_simulation_verysimple.simulation_state
    return model.time_explosion.cgs.value


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
def opacity_state_args(request) -> tuple:
    params = getattr(request, "param", {})
    line_list_nu = np.array(params.get("line_list_nu", [3.999e14, 3.998e14]))
    tau_sobolev = params.get("tau_sobolev", np.zeros((2, 2)))
    no_of_lines = len(line_list_nu)
    no_of_shells = tau_sobolev.shape[1]
    return (
        np.ones(no_of_shells) * 1.0e8,
        np.ones(no_of_shells) * 1.0e4,
        line_list_nu,
        tau_sobolev,
        np.zeros((1, no_of_shells)),
        np.zeros(no_of_lines, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.zeros(1),
        np.zeros((1, no_of_shells)),
        np.zeros(1),
        np.zeros(1),
        np.zeros(1, dtype=np.int64),
        np.zeros((1, no_of_shells)),
        np.zeros(1),
        np.zeros(1),
        np.zeros(no_of_shells),
        np.zeros((1, no_of_shells)),
        np.zeros(1, dtype=np.int64),
        0,
    )


@pytest.fixture
def montecarlo_configuration() -> MonteCarloConfiguration:
    config = MonteCarloConfiguration()
    config.LINE_INTERACTION_TYPE = 0
    config.SURVIVAL_PROBABILITY = 0.0
    config.VPACKET_TAU_RUSSIAN = 10.0
    return config


@pytest.fixture
def characterization_packet(static_packet: RPacket, request) -> RPacket:
    params = getattr(request, "param", {})
    static_packet.nu = 4.0e14
    static_packet.current_shell_id = params.get("current_shell_id", 0)
    static_packet.next_line_id = params.get("next_line_id", 0)
    static_packet.prev_line_id = params.get("prev_line_id", 0)
    return static_packet


@pytest.fixture
def radial_geometry(request) -> NumbaRadial1DGeometry:
    r_outer_first_shell = getattr(request, "param", 8.0e14)
    return NumbaRadial1DGeometry(
        np.array([7.0e14, 8.0e14]),
        np.array([r_outer_first_shell, 3.0e16]),
        np.array([-1.0, -1.0]),
        np.array([-1.0, -1.0]),
    )


@pytest.fixture
def nonhomologous_geometry(request) -> NumbaNonhomologousRadial1DGeometry:
    params = getattr(request, "param", {})
    if params.get("negative_velocity_gradient", False):
        v_inner = np.array([1.5e9, 2.0e9])
        v_outer = np.array([1.0e9, 1.5e9])
    else:
        v_inner = np.array([1.0e9, 1.5e9])
        v_outer = np.array([1.5e9, 2.0e9])

    return NumbaNonhomologousRadial1DGeometry(
        np.array([7.0e14, 8.0e14]),
        np.array([params.get("r_outer_first_shell", 8.0e14), 3.0e16]),
        v_inner,
        v_outer,
    )


@pytest.fixture
def classic_opacity_state(opacity_state_args: tuple) -> OpacityStateNumba:
    return OpacityStateNumba(*opacity_state_args)


@pytest.fixture
def iip_opacity_state(opacity_state_args: tuple) -> OpacityStateNumbaIIP:
    no_of_shells = opacity_state_args[3].shape[1]
    return OpacityStateNumbaIIP(
        *opacity_state_args,
        np.ones((no_of_shells, 1, 1)),
    )


@pytest.fixture
def bulk_estimators():
    return init_estimators_bulk(2)


@pytest.fixture
def line_estimators() -> EstimatorsLine:
    return EstimatorsLine(np.zeros((2, 2)), np.zeros((2, 2)))


@pytest.fixture
def continuum_estimators():
    return init_estimators_continuum((1, 2), 2)


@pytest.fixture
def vpacket_collection() -> VPacketCollection:
    return VPacketCollection(
        source_rpacket_index=0,
        spectrum_frequency_grid=np.array([1.0e14, 2.0e14]),
        number_of_vpackets=0,
        v_packet_spawn_start_frequency=0.0,
        v_packet_spawn_end_frequency=np.inf,
        temporary_v_packet_bins=0,
    )
