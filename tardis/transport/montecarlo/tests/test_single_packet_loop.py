from types import ModuleType

import numpy as np
import numpy.testing as npt
import pytest

import tardis.transport.montecarlo.modes.iip.packet_propagation as iip_propagation
from tardis.conftest import sync_ndarray_assert_allclose
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.model.geometry.radial1d_homologous import (
    NumbaHomologousRadial1DGeometry,
)
from tardis.transport.montecarlo import RPacket
from tardis.transport.montecarlo.configuration.constants import C_SPEED_OF_LIGHT
from tardis.transport.montecarlo.modes.classic import (
    packet_propagation as classic_propagation,
)
from tardis.transport.montecarlo.modes.nonhomologous import (
    packet_propagation as nonhomologous_propagation,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
)

RTOL = 1.0e-12


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

    def __call__(
        self,
        module: ModuleType,
        interactions: list[tuple[float, InteractionType, int]],
    ) -> None:
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
    python_numba_disabled: None,
    patch_common_classic_hooks,
    parametrized_packet: RPacket,
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
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
    patch_common_classic_hooks(
        classic_propagation,
        interactions,
    )

    classic_propagation.packet_propagation(
        parametrized_packet,
        homologous_radial_1d_geometry,
        classic_opacity_state,
        bulk_estimators,
        line_estimators,
        vpacket_collection,
        recording_tracker,
        montecarlo_configuration,
        False,
    )

    assert parametrized_packet.status == PacketStatus.EMITTED
    assert expected_interaction in recording_tracker.events
    assert recording_tracker.events[-1] == InteractionType.BOUNDARY
    sync_ndarray_assert_allclose(
        regression_data,
        np.array(
            [
                parametrized_packet.r,
                parametrized_packet.mu,
                parametrized_packet.nu,
                parametrized_packet.energy,
            ]
        ),
        bulk_estimators.mean_intensity_total,
        bulk_estimators.mean_frequency,
        line_estimators.mean_intensity_blueward,
        line_estimators.energy_deposition_line_rate,
        rtol=RTOL,
    )


def test_classic_packet_propagation_full_relativity_numba_disabled(
    python_numba_disabled: None,
    patch_common_classic_hooks,
    parametrized_packet: RPacket,
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
    classic_opacity_state,
    bulk_estimators,
    line_estimators,
    vpacket_collection,
    recording_tracker: _RecordingTracker,
    montecarlo_configuration,
) -> None:
    patch_common_classic_hooks(
        classic_propagation,
        [
            (0.0, InteractionType.BOUNDARY, 1),
            (0.0, InteractionType.BOUNDARY, 1),
        ],
    )
    initial_mu = parametrized_packet.mu
    initial_nu = parametrized_packet.nu
    initial_energy = parametrized_packet.energy
    velocity = homologous_radial_1d_geometry.get_velocity(
        parametrized_packet.r,
        parametrized_packet.current_shell_id,
    )
    beta = velocity / C_SPEED_OF_LIGHT
    inverse_doppler_factor = (1.0 + initial_mu * beta) / np.sqrt(1.0 - beta**2)
    expected_mu = (initial_mu + beta) / (1.0 + beta * initial_mu)

    classic_propagation.packet_propagation(
        parametrized_packet,
        homologous_radial_1d_geometry,
        classic_opacity_state,
        bulk_estimators,
        line_estimators,
        vpacket_collection,
        recording_tracker,
        montecarlo_configuration,
        True,
    )

    assert parametrized_packet.status == PacketStatus.EMITTED
    npt.assert_allclose(
        parametrized_packet.nu,
        initial_nu * inverse_doppler_factor,
    )
    npt.assert_allclose(
        parametrized_packet.energy,
        initial_energy * inverse_doppler_factor,
    )
    npt.assert_allclose(parametrized_packet.mu, expected_mu)


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
    monkeypatch: pytest.MonkeyPatch,
    patch_common_classic_hooks,
    parametrized_packet: RPacket,
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
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
    patch_common_classic_hooks(
        iip_propagation,
        interactions,
    )
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
        homologous_radial_1d_geometry,
        iip_opacity_state,
        bulk_estimators,
        line_estimators,
        continuum_estimators,
        recording_tracker,
        montecarlo_configuration,
        True,
    )

    assert parametrized_packet.status == PacketStatus.EMITTED
    assert first_interaction in recording_tracker.events
    assert recording_tracker.events[-1] == InteractionType.BOUNDARY
    sync_ndarray_assert_allclose(
        regression_data,
        np.array(
            [
                parametrized_packet.r,
                parametrized_packet.mu,
                parametrized_packet.nu,
                parametrized_packet.energy,
            ]
        ),
        bulk_estimators.mean_intensity_total,
        bulk_estimators.mean_frequency,
        line_estimators.mean_intensity_blueward,
        line_estimators.energy_deposition_line_rate,
        continuum_estimators.photo_ion_estimator,
        rtol=RTOL,
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
    python_numba_disabled: None,
    patch_common_classic_hooks,
    parametrized_packet: RPacket,
    radial_1d_geometry: NumbaRadial1DGeometry,
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
    patch_common_classic_hooks(
        nonhomologous_propagation,
        interactions,
    )

    nonhomologous_propagation.packet_propagation(
        parametrized_packet,
        radial_1d_geometry,
        classic_opacity_state,
        bulk_estimators,
        line_estimators,
        vpacket_collection,
        recording_tracker,
        montecarlo_configuration,
        False,
    )

    assert parametrized_packet.status == PacketStatus.EMITTED
    assert first_interaction in recording_tracker.events
    assert recording_tracker.events[-1] == InteractionType.BOUNDARY
    sync_ndarray_assert_allclose(
        regression_data,
        np.array(
            [
                parametrized_packet.r,
                parametrized_packet.mu,
                parametrized_packet.nu,
                parametrized_packet.energy,
            ]
        ),
        bulk_estimators.mean_intensity_total,
        bulk_estimators.mean_frequency,
        line_estimators.mean_intensity_blueward,
        line_estimators.energy_deposition_line_rate,
        rtol=RTOL,
    )
