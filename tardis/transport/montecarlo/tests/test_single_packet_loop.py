import os

import numpy as np
import numpy.testing as npt
import pytest

import tardis.transport.montecarlo.modes.classic.packet_propagation as classic_propagation
import tardis.transport.montecarlo.modes.iip.packet_propagation as iip_propagation
import tardis.transport.montecarlo.modes.nonhomologous.packet_propagation as nonhomologous_propagation
from tardis.conftest import assert_synced_allclose
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.model.geometry.radial1d_nonhomologous import (
    NumbaNonhomologousRadial1DGeometry,
)
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.opacities.opacity_state_numba_iip import OpacityStateNumbaIIP
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
from tardis.transport.montecarlo.modes.classic.packet_propagation import (
    packet_propagation,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    VPacketCollection,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
)
from tardis.transport.montecarlo.tests.test_transport_characterization import (
    _opacity_state_args,
)


class _RecordingTracker:
    def __init__(self) -> None:
        self.events: list[InteractionType] = []
        self.boundaries: list[tuple[int, int]] = []

    def track_boundary_event(
        self, r_packet, from_shell_id: int, to_shell_id: int
    ) -> None:
        self.events.append(InteractionType.BOUNDARY)
        self.boundaries.append((from_shell_id, to_shell_id))

    def track_line_interaction_before(self, r_packet) -> None:
        self.events.append(InteractionType.LINE)

    def track_line_interaction_after(self, r_packet) -> None:
        pass

    def track_escattering_interaction_before(self, r_packet) -> None:
        self.events.append(InteractionType.ESCATTERING)

    def track_escattering_interaction_after(self, r_packet) -> None:
        pass

    def track_continuum_interaction_before(self, r_packet) -> None:
        self.events.append(InteractionType.CONTINUUM_PROCESS)

    def track_continuum_interaction_after(self, r_packet) -> None:
        pass


def _skip_unless_python_propagation() -> None:
    if os.environ.get("NUMBA_DISABLE_JIT") != "1":
        pytest.skip(
            "packet_propagation dispatch characterization uses monkeypatching "
            "and only runs with NUMBA_DISABLE_JIT=1"
        )


@pytest.fixture
def montecarlo_configuration() -> MonteCarloConfiguration:
    config = MonteCarloConfiguration()
    config.LINE_INTERACTION_TYPE = 0
    config.SURVIVAL_PROBABILITY = 0.0
    config.VPACKET_TAU_RUSSIAN = 10.0
    return config


@pytest.fixture
def characterization_packet() -> RPacket:
    packet = RPacket(
        r=7.5e14,
        mu=0.3,
        nu=4.0e14,
        energy=0.9,
        seed=1963,
        index=0,
    )
    packet.current_shell_id = 0
    packet.next_line_id = 0
    packet.prev_line_id = 0
    return packet


@pytest.fixture
def radial_geometry() -> NumbaRadial1DGeometry:
    return NumbaRadial1DGeometry(
        np.array([7.0e14, 8.0e14]),
        np.array([8.0e14, 3.0e16]),
        np.array([-1.0, -1.0]),
        np.array([-1.0, -1.0]),
    )


@pytest.fixture
def nonhomologous_geometry() -> NumbaNonhomologousRadial1DGeometry:
    return NumbaNonhomologousRadial1DGeometry(
        np.array([7.0e14, 8.0e14]),
        np.array([8.0e14, 3.0e16]),
        np.array([1.0e9, 1.5e9]),
        np.array([1.5e9, 2.0e9]),
    )


@pytest.fixture
def classic_opacity_state() -> OpacityStateNumba:
    return OpacityStateNumba(
        *_opacity_state_args(np.array([3.999e14, 3.998e14]), np.zeros((2, 2)))
    )


@pytest.fixture
def iip_opacity_state() -> OpacityStateNumbaIIP:
    return OpacityStateNumbaIIP(
        *_opacity_state_args(np.array([3.999e14, 3.998e14]), np.zeros((2, 2))),
        np.ones((2, 1, 1)),
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
def recording_tracker() -> _RecordingTracker:
    return _RecordingTracker()


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


def _patch_common_classic_hooks(monkeypatch, module, interactions) -> None:
    interaction_queue = list(interactions)

    def fake_trace_packet(*args, **kwargs):
        assert interaction_queue
        return interaction_queue.pop(0)

    monkeypatch.setattr(module, "trace_packet", fake_trace_packet)
    monkeypatch.setattr(
        module, "trace_vpacket_volley", lambda *args: None, raising=False
    )
    monkeypatch.setattr(
        module, "chi_electron_calculator", lambda *args: 1.0e-20
    )
    monkeypatch.setattr(
        module,
        "line_scatter_event",
        lambda packet, *args, **kwargs: setattr(packet, "next_line_id", 1),
    )
    monkeypatch.setattr(
        module,
        "thomson_scatter",
        lambda packet, *args, **kwargs: setattr(packet, "mu", -packet.mu),
    )


def _assert_dispatch_result(
    packet: RPacket,
    tracker: _RecordingTracker,
    expected_interaction: InteractionType,
) -> None:
    assert packet.status == PacketStatus.EMITTED
    assert expected_interaction in tracker.events
    assert tracker.events[-1] == InteractionType.BOUNDARY


@pytest.mark.parametrize(
    ("first_interaction", "expected_interaction"),
    [
        (InteractionType.BOUNDARY, InteractionType.BOUNDARY),
        (InteractionType.LINE, InteractionType.LINE),
        (InteractionType.ESCATTERING, InteractionType.ESCATTERING),
    ],
)
def test_classic_packet_propagation_dispatch_characterization(
    monkeypatch,
    characterization_packet: RPacket,
    radial_geometry,
    classic_opacity_state,
    bulk_estimators,
    line_estimators,
    vpacket_collection: VPacketCollection,
    recording_tracker: _RecordingTracker,
    montecarlo_configuration: MonteCarloConfiguration,
    first_interaction: InteractionType,
    expected_interaction: InteractionType,
    regression_data,
) -> None:
    _skip_unless_python_propagation()
    interactions = [(1.0e12, first_interaction, 1)]
    if first_interaction != InteractionType.BOUNDARY:
        interactions.append((1.0e12, InteractionType.BOUNDARY, 1))
    interactions.append((1.0e12, InteractionType.BOUNDARY, 1))
    _patch_common_classic_hooks(monkeypatch, classic_propagation, interactions)

    classic_propagation.packet_propagation(
        characterization_packet,
        radial_geometry,
        5.2e7,
        classic_opacity_state,
        bulk_estimators,
        line_estimators,
        vpacket_collection,
        recording_tracker,
        montecarlo_configuration,
    )

    _assert_dispatch_result(
        characterization_packet,
        recording_tracker,
        expected_interaction,
    )
    assert_synced_allclose(
        regression_data,
        np.array(
            [
                characterization_packet.r,
                characterization_packet.mu,
                characterization_packet.nu,
                characterization_packet.energy,
            ]
        ),
        bulk_estimators.mean_intensity_total,
        bulk_estimators.mean_frequency,
        line_estimators.mean_intensity_blueward,
        line_estimators.energy_deposition_line_rate,
    )


@pytest.mark.parametrize(
    "first_interaction",
    [
        InteractionType.LINE,
        InteractionType.ESCATTERING,
        InteractionType.CONTINUUM_PROCESS,
    ],
)
def test_iip_packet_propagation_dispatch_characterization(
    monkeypatch,
    characterization_packet: RPacket,
    radial_geometry,
    iip_opacity_state,
    bulk_estimators,
    line_estimators,
    continuum_estimators,
    recording_tracker: _RecordingTracker,
    montecarlo_configuration: MonteCarloConfiguration,
    first_interaction: InteractionType,
    regression_data,
) -> None:
    _skip_unless_python_propagation()
    interactions = [
        (1.0e12, first_interaction, 1),
        (1.0e12, InteractionType.BOUNDARY, 1),
        (1.0e12, InteractionType.BOUNDARY, 1),
    ]
    _patch_common_classic_hooks(monkeypatch, iip_propagation, interactions)
    monkeypatch.setattr(
        iip_propagation,
        "chi_continuum_calculator",
        lambda *args: (
            1.0,
            np.array([1.0]),
            np.array([0], dtype=np.int64),
            np.array([1.0]),
            0.0,
        ),
    )
    monkeypatch.setattr(
        iip_propagation, "update_estimators_bound_free", lambda *args: None
    )
    monkeypatch.setattr(
        iip_propagation, "continuum_event", lambda *args, **kwargs: None
    )

    iip_propagation.packet_propagation(
        characterization_packet,
        radial_geometry,
        5.2e7,
        iip_opacity_state,
        bulk_estimators,
        line_estimators,
        continuum_estimators,
        recording_tracker,
        montecarlo_configuration,
    )

    _assert_dispatch_result(
        characterization_packet,
        recording_tracker,
        first_interaction,
    )
    assert_synced_allclose(
        regression_data,
        np.array(
            [
                characterization_packet.r,
                characterization_packet.mu,
                characterization_packet.nu,
                characterization_packet.energy,
            ]
        ),
        bulk_estimators.mean_intensity_total,
        bulk_estimators.mean_frequency,
        line_estimators.mean_intensity_blueward,
        line_estimators.energy_deposition_line_rate,
        continuum_estimators.photo_ion_estimator,
    )


@pytest.mark.parametrize(
    "first_interaction",
    [
        InteractionType.BOUNDARY,
        InteractionType.LINE,
        InteractionType.ESCATTERING,
    ],
)
def test_nonhomologous_packet_propagation_dispatch_characterization(
    monkeypatch,
    characterization_packet: RPacket,
    nonhomologous_geometry,
    classic_opacity_state,
    bulk_estimators,
    line_estimators,
    vpacket_collection: VPacketCollection,
    recording_tracker: _RecordingTracker,
    montecarlo_configuration: MonteCarloConfiguration,
    first_interaction: InteractionType,
    regression_data,
) -> None:
    _skip_unless_python_propagation()
    interactions = [(1.0e12, first_interaction, 1)]
    if first_interaction != InteractionType.BOUNDARY:
        interactions.append((1.0e12, InteractionType.BOUNDARY, 1))
    interactions.append((1.0e12, InteractionType.BOUNDARY, 1))
    _patch_common_classic_hooks(
        monkeypatch, nonhomologous_propagation, interactions
    )

    nonhomologous_propagation.packet_propagation(
        characterization_packet,
        nonhomologous_geometry,
        classic_opacity_state,
        bulk_estimators,
        line_estimators,
        vpacket_collection,
        recording_tracker,
        montecarlo_configuration,
    )

    _assert_dispatch_result(
        characterization_packet,
        recording_tracker,
        first_interaction,
    )
    assert_synced_allclose(
        regression_data,
        np.array(
            [
                characterization_packet.r,
                characterization_packet.mu,
                characterization_packet.nu,
                characterization_packet.energy,
            ]
        ),
        bulk_estimators.mean_intensity_total,
        bulk_estimators.mean_frequency,
        line_estimators.mean_intensity_blueward,
        line_estimators.energy_deposition_line_rate,
    )


@pytest.mark.xfail(
    reason="Need to update for new mode architecture with correct parameters"
)
# TODO: Update test to provide all required parameters for packet_propagation
def test_verysimple_single_packet_loop(
    verysimple_numba_radial_1d_geometry,
    verysimple_time_explosion,
    verysimple_opacity_state,
    verysimple_estimators,
    verysimple_vpacket_collection,
    verysimple_packet_collection,
):
    pytest.skip("Test needs to be updated for new mode architecture")
    numba_radial_1d_geometry = verysimple_numba_radial_1d_geometry
    packet_collection = verysimple_packet_collection
    vpacket_collection = verysimple_vpacket_collection
    time_explosion = verysimple_time_explosion
    opacity_state = verysimple_opacity_state
    numba_estimators = verysimple_estimators

    i = 0
    r_packet = RPacket(
        numba_radial_1d_geometry.r_inner[0],
        packet_collection.packets_input_mu[i],
        packet_collection.packets_input_nu[i],
        packet_collection.packets_input_energy[i],
        i,
    )
    # packet_propagation requires: r_packet, geometry, time_explosion, opacity_state,
    # estimators_bulk, estimators_line, vpacket_collection, rpacket_tracker, montecarlo_configuration
    # This test needs to be updated with all required parameters
    packet_propagation(
        r_packet,
        numba_radial_1d_geometry,
        time_explosion,
        opacity_state,
        numba_estimators,
        vpacket_collection,
    )

    npt.assert_almost_equal(r_packet.nu, 1053057938883272.8)
    npt.assert_almost_equal(r_packet.mu, 0.9611146425440562)
    npt.assert_almost_equal(r_packet.energy, 0.10327717505563379)


@pytest.mark.xfail(reason="To be implemented")
def test_set_packet_props_partial_relativity():
    raise AssertionError()


@pytest.mark.xfail(reason="To be implemented")
def test_set_packet_props_full_relativity():
    raise AssertionError()
