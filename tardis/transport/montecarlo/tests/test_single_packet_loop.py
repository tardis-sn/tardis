import numpy as np
import numpy.testing as npt
import pytest

import tardis.transport.montecarlo.modes.iip.packet_propagation as iip_propagation
from tardis.conftest import assert_synced_allclose
from tardis.transport.montecarlo import RPacket
from tardis.transport.montecarlo.modes.classic import (
    packet_propagation as classic_propagation,
)
from tardis.transport.montecarlo.modes.classic.packet_propagation import (
    packet_propagation,
)
from tardis.transport.montecarlo.modes.nonhomologous import (
    packet_propagation as nonhomologous_propagation,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
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


class _CommonPacketPropagationPatcher:
    def __init__(self, monkeypatch) -> None:
        self.monkeypatch = monkeypatch
        self.interaction_queue = []

    def trace_packet(self, *args, **kwargs):
        assert self.interaction_queue
        return self.interaction_queue.pop(0)

    def __call__(self, module, interactions) -> None:
        self.interaction_queue = list(interactions)
        self.monkeypatch.setattr(module, "trace_packet", self.trace_packet)
        self.monkeypatch.setattr(
            module, "trace_vpacket_volley", lambda *args: None, raising=False
        )
        self.monkeypatch.setattr(
            module, "chi_electron_calculator", lambda *args: 1.0e-20
        )
        self.monkeypatch.setattr(
            module,
            "line_scatter_event",
            lambda packet, *args, **kwargs: setattr(packet, "next_line_id", 1),
        )
        self.monkeypatch.setattr(
            module,
            "thomson_scatter",
            lambda packet, *args, **kwargs: setattr(packet, "mu", -packet.mu),
        )


@pytest.fixture
def recording_tracker() -> _RecordingTracker:
    return _RecordingTracker()


@pytest.fixture
def patch_common_classic_hooks(monkeypatch):
    return _CommonPacketPropagationPatcher(monkeypatch)


@pytest.mark.parametrize(
    ("first_interaction", "expected_interaction"),
    [
        (InteractionType.BOUNDARY, InteractionType.BOUNDARY),
        (InteractionType.LINE, InteractionType.LINE),
        (InteractionType.ESCATTERING, InteractionType.ESCATTERING),
    ],
)
def test_classic_packet_propagation_dispatch_numba_disabled(
    python_numba_disabled,
    patch_common_classic_hooks,
    parametrized_packet: RPacket,
    radial_geometry,
    classic_opacity_state,
    bulk_estimators,
    line_estimators,
    vpacket_collection,
    recording_tracker: _RecordingTracker,
    montecarlo_configuration,
    first_interaction: InteractionType,
    expected_interaction: InteractionType,
    regression_data,
) -> None:
    # Force the propagation loop through one selected interaction, then an
    # outer boundary, so this characterizes dispatch bookkeeping directly.
    interactions = [(1.0e12, first_interaction, 1)]
    if first_interaction != InteractionType.BOUNDARY:
        interactions.append((1.0e12, InteractionType.BOUNDARY, 1))
    interactions.append((1.0e12, InteractionType.BOUNDARY, 1))
    patch_common_classic_hooks(classic_propagation, interactions)

    classic_propagation.packet_propagation(
        parametrized_packet,
        radial_geometry,
        5.2e7,
        classic_opacity_state,
        bulk_estimators,
        line_estimators,
        vpacket_collection,
        recording_tracker,
        montecarlo_configuration,
    )

    assert parametrized_packet.status == PacketStatus.EMITTED
    assert expected_interaction in recording_tracker.events
    assert recording_tracker.events[-1] == InteractionType.BOUNDARY
    assert_synced_allclose(
        regression_data,
        np.array([
            parametrized_packet.r,
            parametrized_packet.mu,
            parametrized_packet.nu,
            parametrized_packet.energy,
        ]),
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
def test_iip_packet_propagation_dispatch_numba_disabled(
    python_numba_disabled,
    monkeypatch,
    patch_common_classic_hooks,
    parametrized_packet: RPacket,
    radial_geometry,
    iip_opacity_state,
    bulk_estimators,
    line_estimators,
    continuum_estimators,
    recording_tracker: _RecordingTracker,
    montecarlo_configuration,
    first_interaction: InteractionType,
    regression_data,
) -> None:
    # Force IIP dispatch paths independently of continuum opacity details.
    interactions = [
        (1.0e12, first_interaction, 1),
        (1.0e12, InteractionType.BOUNDARY, 1),
        (1.0e12, InteractionType.BOUNDARY, 1),
    ]
    patch_common_classic_hooks(iip_propagation, interactions)
    # Continuum handling has extra collaborators; patch them only enough to
    # make the dispatch branch observable.
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
        parametrized_packet,
        radial_geometry,
        5.2e7,
        iip_opacity_state,
        bulk_estimators,
        line_estimators,
        continuum_estimators,
        recording_tracker,
        montecarlo_configuration,
    )

    assert parametrized_packet.status == PacketStatus.EMITTED
    assert first_interaction in recording_tracker.events
    assert recording_tracker.events[-1] == InteractionType.BOUNDARY
    assert_synced_allclose(
        regression_data,
        np.array([
            parametrized_packet.r,
            parametrized_packet.mu,
            parametrized_packet.nu,
            parametrized_packet.energy,
        ]),
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
def test_nonhomologous_packet_propagation_dispatch_numba_disabled(
    python_numba_disabled,
    patch_common_classic_hooks,
    parametrized_packet: RPacket,
    nonhomologous_geometry,
    classic_opacity_state,
    bulk_estimators,
    line_estimators,
    vpacket_collection,
    recording_tracker: _RecordingTracker,
    montecarlo_configuration,
    first_interaction: InteractionType,
    regression_data,
) -> None:
    # Force the nonhomologous loop through one selected interaction, then an
    # outer boundary, so this characterizes dispatch bookkeeping directly.
    interactions = [(1.0e12, first_interaction, 1)]
    if first_interaction != InteractionType.BOUNDARY:
        interactions.append((1.0e12, InteractionType.BOUNDARY, 1))
    interactions.append((1.0e12, InteractionType.BOUNDARY, 1))
    patch_common_classic_hooks(nonhomologous_propagation, interactions)

    nonhomologous_propagation.packet_propagation(
        parametrized_packet,
        nonhomologous_geometry,
        classic_opacity_state,
        bulk_estimators,
        line_estimators,
        vpacket_collection,
        recording_tracker,
        montecarlo_configuration,
    )

    assert parametrized_packet.status == PacketStatus.EMITTED
    assert first_interaction in recording_tracker.events
    assert recording_tracker.events[-1] == InteractionType.BOUNDARY
    assert_synced_allclose(
        regression_data,
        np.array([
            parametrized_packet.r,
            parametrized_packet.mu,
            parametrized_packet.nu,
            parametrized_packet.energy,
        ]),
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
    # estimators_bulk, estimators_line, vpacket_collection, rpacket_tracker,
    # montecarlo_configuration
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
