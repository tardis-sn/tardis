from collections.abc import Callable
from types import SimpleNamespace

import numpy as np
import numpy.testing as npt
import pytest

import tardis.transport.montecarlo.modes.iip.packet_propagation as iip_propagation
from tardis import constants as const
from tardis.conftest import sync_ndarray_assert_allclose
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.model.geometry.radial1d_homologous import (
    NumbaHomologousRadial1DGeometry,
)
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.montecarlo import RPacket
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.configuration.constants import SIGMA_THOMSON
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
)
from tardis.transport.montecarlo.modes.classic import (
    packet_propagation as classic_propagation,
)
from tardis.transport.montecarlo.modes.classic.packet_propagation import (
    packet_propagation,
)
from tardis.transport.montecarlo.modes.nonhomologous import (
    interaction_event_callers as nonhomologous_callers,
)
from tardis.transport.montecarlo.modes.nonhomologous import (
    packet_propagation as nonhomologous_propagation,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full import (
    TrackerFull,
)

RTOL = 1.0e-12
C_CGS = const.c.cgs.value


@pytest.mark.parametrize(
    "mode", ["homologous_partial", "homologous_full", "nonhomologous_partial"]
)
def test_transparent_ejecta_escapes_at_geometric_boundary(
    mode: str,
    parametrized_packet: RPacket,
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
    classic_opacity_state: OpacityStateNumba,
    bulk_estimators: EstimatorsBulk,
    line_estimators: EstimatorsLine,
    montecarlo_configuration: MonteCarloConfiguration,
    set_seed_fixture: Callable[[int], None],
) -> None:
    """Verify transparent transport against a spherical chord (verification).

    Claim: zero-opacity packets escape at the outer sphere, accumulating only
    the per-segment bulk path estimators. Regime: two shells, outward oblique
    ray, zero electron density/tau, retained lines outside the frequency sweep.
    Verification: impact-parameter geometry and analytic Lorentz factors,
    independently of the movement/frame/estimator helpers. The bulk oracle
    characterizes the current segment-start Doppler quadrature, rather than
    claiming an exact continuous integral. rtol=1e-12, atol=0 allows floating
    point geometry/transform rounding. Detects wrong frame, movement, shell
    dispatch, or residual line/electron opacity.
    """
    geometry = homologous_radial_1d_geometry
    time_explosion = geometry.time_explosion
    full = mode == "homologous_full"
    propagate = classic_propagation.packet_propagation
    if mode == "nonhomologous_partial":
        geometry = NumbaRadial1DGeometry(
            geometry.r_inner,
            geometry.r_outer,
            np.array([1.0e8, 2.0e8]),
            np.array([2.0e8, 3.0e8]),
        )
        propagate = nonhomologous_propagation.packet_propagation
    montecarlo_configuration.ENABLE_FULL_RELATIVITY = full
    classic_opacity_state.electron_density[:] = 0.0
    classic_opacity_state.tau_sobolev[:] = 0.0
    # Both retained lines lie blueward of the entire escaping packet sweep.
    classic_opacity_state.line_list_nu[:] = [1.0e16, 9.0e15]
    packet = parametrized_packet
    radius, mu, energy, frequency = (
        packet.r,
        packet.mu,
        packet.energy,
        packet.nu,
    )
    if mode == "nonhomologous_partial":
        velocity = 1.0e8 + 1.0e8 * (radius - geometry.r_inner[0]) / (
            geometry.r_outer[0] - geometry.r_inner[0]
        )
    else:
        velocity = radius / time_explosion
    beta = velocity / C_CGS
    inverse_doppler = 1.0 / (1.0 - beta * mu)
    if full:
        inverse_doppler = (1.0 + beta * mu) / np.sqrt(1.0 - beta**2)
        mu = (mu + beta) / (1.0 + beta * mu)
    lab_energy, lab_frequency = (
        energy * inverse_doppler,
        frequency * inverse_doppler,
    )
    impact_squared = radius**2 * (1.0 - mu**2)
    endpoints = np.sqrt(geometry.r_outer**2 - impact_squared)
    starts = np.array([radius * mu, endpoints[0]])
    lengths = endpoints - starts
    start_radii = np.array([radius, geometry.r_outer[0]])
    start_mus = starts / start_radii
    velocities = start_radii / time_explosion
    if mode == "nonhomologous_partial":
        velocities = np.array([velocity, 2.0e8])
    betas = velocities / C_CGS
    dopplers = 1.0 - betas * start_mus
    if full:
        dopplers /= np.sqrt(1.0 - betas**2)
    comoving_lengths = lengths * (dopplers if full else 1.0)
    expected_intensity = lab_energy * dopplers * comoving_lengths
    expected_frequency = expected_intensity * lab_frequency * dopplers
    tracker = TrackerFull(8, 2)
    set_seed_fixture(1963)
    propagate(
        packet,
        geometry,
        time_explosion,
        classic_opacity_state,
        bulk_estimators,
        line_estimators,
        tracker,
        montecarlo_configuration,
    )
    tracker.finalize()

    npt.assert_allclose(packet.r, geometry.r_outer[-1], rtol=RTOL, atol=0)
    npt.assert_allclose(
        packet.mu, endpoints[-1] / geometry.r_outer[-1], rtol=RTOL, atol=0
    )
    npt.assert_allclose(
        [packet.energy, packet.nu],
        [lab_energy, lab_frequency],
        rtol=RTOL,
        atol=0,
    )
    assert packet.current_shell_id == 1
    assert packet.status == PacketStatus.EMITTED
    npt.assert_array_equal(
        tracker.interaction_type, np.full(4, InteractionType.BOUNDARY)
    )
    npt.assert_array_equal(tracker.before_shell_id, [-1, 0, 1, 1])
    npt.assert_array_equal(tracker.after_shell_id, [0, 1, 2, 2])
    npt.assert_allclose(
        tracker.radius,
        [radius, *geometry.r_outer, geometry.r_outer[-1]],
        rtol=RTOL,
        atol=0,
    )
    npt.assert_allclose(
        bulk_estimators.mean_intensity_total,
        expected_intensity,
        rtol=RTOL,
        atol=0,
    )
    npt.assert_allclose(
        bulk_estimators.mean_frequency, expected_frequency, rtol=RTOL, atol=0
    )
    npt.assert_array_equal(line_estimators.mean_intensity_blueward, 0.0)
    npt.assert_array_equal(line_estimators.energy_deposition_line_rate, 0.0)


@pytest.mark.parametrize(
    "mode", ["homologous_partial", "homologous_full", "nonhomologous_partial"]
)
@pytest.mark.parametrize(
    ("tau", "expected_line_events"), [(0.0, 0), (100.0, 1)]
)
def test_isolated_line_obeys_zero_and_opaque_limits(
    mode: str,
    tau: float,
    expected_line_events: int,
    parametrized_packet: RPacket,
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
    classic_opacity_state: OpacityStateNumba,
    bulk_estimators: EstimatorsBulk,
    line_estimators: EstimatorsLine,
    montecarlo_configuration: MonteCarloConfiguration,
    set_seed_fixture: Callable[[int], None],
) -> None:
    """A tau=0 line is never absorbed; a tau=100 line is absorbed exactly once.

    Regime: electron-free two-shell ejecta, one reachable line plus a far
    red sentinel, SCATTER. The interaction probability is 1 - exp(-tau): 0 and
    1 - 4e-44, so the outcome does not depend on the random seed.
    """
    geometry = homologous_radial_1d_geometry
    time_explosion = geometry.time_explosion
    propagate = classic_propagation.packet_propagation
    if mode == "nonhomologous_partial":
        geometry = NumbaRadial1DGeometry(
            geometry.r_inner,
            geometry.r_outer,
            geometry.v_inner,
            geometry.v_outer,
        )
        propagate = nonhomologous_propagation.packet_propagation
    montecarlo_configuration.ENABLE_FULL_RELATIVITY = mode == "homologous_full"
    opacity = classic_opacity_state
    opacity.electron_density[:] = 0.0
    opacity.line_list_nu[:] = [3.9999e14, 1.0e9]
    opacity.tau_sobolev[:] = 0.0
    opacity.tau_sobolev[0, :] = tau
    tracker = TrackerFull(16, 2)
    set_seed_fixture(1963)
    propagate(
        parametrized_packet,
        geometry,
        time_explosion,
        opacity,
        bulk_estimators,
        line_estimators,
        tracker,
        montecarlo_configuration,
    )
    tracker.finalize()

    line_events = tracker.interaction_type == InteractionType.LINE
    assert np.count_nonzero(line_events) == expected_line_events
    npt.assert_array_equal(tracker.line_absorb_id[line_events], 0)
    npt.assert_array_equal(tracker.line_emit_id[line_events], 0)


@pytest.mark.parametrize(
    "mode", ["homologous_partial", "homologous_full", "nonhomologous_partial"]
)
def test_electron_opaque_limit_selects_thomson_scattering(
    mode: str,
    parametrized_packet: RPacket,
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
    classic_opacity_state: OpacityStateNumba,
    bulk_estimators: EstimatorsBulk,
    line_estimators: EstimatorsLine,
    montecarlo_configuration: MonteCarloConfiguration,
    set_seed_fixture: Callable[[int], None],
) -> None:
    """Verify electron opacity selects the exponential free path (verification).

    Claim: the first physical event in an electron-opaque shell is Thomson
    scattering at s=-log(U)/chi. Regime: shell optical depth 100, zero line
    tau and retained blue lines outside the sweep. Verification: independent
    spherical chord and 1-exp(-100) probability, including the FULL Doppler
    opacity factor. Exact event types and rtol=1e-12, atol=0 for geometry
    detect missing electron opacity or an incorrectly wired process/frame.
    """
    geometry = homologous_radial_1d_geometry
    time_explosion = geometry.time_explosion
    propagate = classic_propagation.packet_propagation
    if mode == "nonhomologous_partial":
        geometry = NumbaRadial1DGeometry(
            geometry.r_inner,
            geometry.r_outer,
            geometry.v_inner,
            geometry.v_outer,
        )
        propagate = nonhomologous_propagation.packet_propagation
    full = mode == "homologous_full"
    montecarlo_configuration.ENABLE_FULL_RELATIVITY = full
    packet = parametrized_packet
    radius, mu = packet.r, packet.mu
    beta = radius / time_explosion / C_CGS
    if full:
        mu = (mu + beta) / (1.0 + beta * mu)
    chord = (
        np.sqrt(geometry.r_outer[0] ** 2 - radius**2 * (1.0 - mu**2))
        - radius * mu
    )
    doppler = (1.0 - beta * mu) / np.sqrt(1.0 - beta**2) if full else 1.0
    opacity = classic_opacity_state
    opacity.electron_density[:] = 0.0
    opacity.electron_density[0] = 100.0 / (chord * SIGMA_THOMSON * doppler)
    opacity.tau_sobolev[:] = 0.0
    opacity.line_list_nu[:] = [1.0e16, 9.0e15]
    seed = 1963
    uniform = np.random.RandomState(seed).random_sample()
    assert uniform > np.exp(-100.0)
    distance = -np.log(uniform) * chord / 100.0
    expected_radius = np.sqrt(
        radius**2 + distance**2 + 2.0 * radius * distance * mu
    )
    tracker = TrackerFull(64, 2)
    set_seed_fixture(seed)
    propagate(
        packet,
        geometry,
        time_explosion,
        opacity,
        bulk_estimators,
        line_estimators,
        tracker,
        montecarlo_configuration,
    )
    tracker.finalize()
    assert tracker.interaction_type[1] == InteractionType.ESCATTERING
    npt.assert_allclose(tracker.radius[1], expected_radius, rtol=RTOL, atol=0)
    assert not np.any(tracker.interaction_type == InteractionType.LINE)
    npt.assert_array_equal(line_estimators.mean_intensity_blueward, 0.0)
    npt.assert_array_equal(line_estimators.energy_deposition_line_rate, 0.0)


class _RecordingTracker:
    def __init__(self) -> None:
        self.events: list[InteractionType] = []
        self.boundaries: list[tuple[int, int]] = []
        self.max_line_events: int | None = None

    def track_boundary_event(
        self, r_packet, from_shell_id: int, to_shell_id: int
    ) -> None:
        self.events.append(InteractionType.BOUNDARY)
        self.boundaries.append((from_shell_id, to_shell_id))

    def track_line_interaction_before(self, r_packet: RPacket) -> None:
        if self.max_line_events is not None:
            line_event_count = self.events.count(InteractionType.LINE)
            assert line_event_count < self.max_line_events, (
                "Packet repeated a line before leaving the shell"
            )
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
    def __init__(self, monkeypatch: pytest.MonkeyPatch) -> None:
        self.monkeypatch = monkeypatch
        self.interaction_queue = []
        self.trace_inputs: list[tuple[int, int]] = []

    def trace_packet(
        self, r_packet: RPacket, *args: object, **kwargs: object
    ) -> tuple[object, ...]:
        assert self.interaction_queue
        self.trace_inputs.append((r_packet.next_line_id, r_packet.prev_line_id))
        return self.interaction_queue.pop(0)

    def line_scatter_event(
        self, r_packet: RPacket, *args: object, **kwargs: object
    ) -> None:
        # The mocked handler emits ID 0, excluding it in both directions.
        r_packet.next_line_id = 1
        r_packet.prev_line_id = -1

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
            self.line_scatter_event,
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


@pytest.fixture
def emitted_line_cursors(
    patch_common_classic_hooks: _CommonPacketPropagationPatcher,
) -> _CommonPacketPropagationPatcher:
    """Record trace-entry cursors after the shared handler emits line 0."""
    interactions = [
        (1.0e12, InteractionType.LINE, 1),
        (1.0e12, InteractionType.BOUNDARY, 1),
        (1.0e12, InteractionType.BOUNDARY, 1),
    ]
    patch_common_classic_hooks(nonhomologous_propagation, interactions)
    return patch_common_classic_hooks


@pytest.fixture
def scatter_packet_loop(
    monkeypatch: pytest.MonkeyPatch,
    recording_tracker: _RecordingTracker,
    montecarlo_configuration: MonteCarloConfiguration,
) -> Callable[..., None]:
    """Use real line scattering with radial emission and one line event."""
    recording_tracker.max_line_events = 1
    montecarlo_configuration.LINE_INTERACTION_TYPE = 0
    montecarlo_configuration.DISABLE_LINE_SCATTERING = False
    monkeypatch.setattr(nonhomologous_callers, "get_random_mu", lambda: 1.0)
    monkeypatch.setattr(
        nonhomologous_propagation,
        "line_scatter_event",
        getattr(
            nonhomologous_callers.line_scatter_event,
            "py_func",
            nonhomologous_callers.line_scatter_event,
        ),
    )
    return getattr(
        nonhomologous_propagation.packet_propagation,
        "py_func",
        nonhomologous_propagation.packet_propagation,
    )


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
        homologous_radial_1d_geometry,
        5.2e7,
        classic_opacity_state,
        bulk_estimators,
        line_estimators,
        recording_tracker,
        montecarlo_configuration,
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
        homologous_radial_1d_geometry,
        5.2e7,
        iip_opacity_state,
        SimpleNamespace(bf_threshold_list_nu=np.array([])),
        bulk_estimators,
        line_estimators,
        continuum_estimators,
        recording_tracker,
        montecarlo_configuration,
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
    recording_tracker: _RecordingTracker,
    montecarlo_configuration,
    first_interaction: InteractionType,
    regression_data,
) -> None:
    """Claim: Nonhomologous propagation dispatches and records each event type.

    Regime: Numba-disabled propagation with mocked boundary, line, or electron
    events followed by escape; interaction physics is mocked.
    Verification: Check the dispatched event, final boundary, and emitted
    status, then compare packet state and estimators with regression data.
    This provides orchestration and regression evidence.
    """
    # Force the nonhomologous loop through one selected interaction, then an
    # outer boundary, so this characterizes dispatch bookkeeping directly.
    interactions = [(1.0e12, first_interaction, 1)]
    if first_interaction != InteractionType.BOUNDARY:
        interactions.append((1.0e12, InteractionType.BOUNDARY, 1))
    interactions.append((1.0e12, InteractionType.BOUNDARY, 1))
    patch_common_classic_hooks(nonhomologous_propagation, interactions)

    nonhomologous_propagation.packet_propagation(
        parametrized_packet,
        radial_1d_geometry,
        0.0,
        classic_opacity_state,
        bulk_estimators,
        line_estimators,
        recording_tracker,
        montecarlo_configuration,
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


@pytest.mark.parametrize(
    "radial_1d_geometry", [{"negative_velocity_gradient": True}], indirect=True
)
def test_nonhomologous_packet_loop_preserves_emitted_line_cursors(
    python_numba_disabled: None,
    emitted_line_cursors: _CommonPacketPropagationPatcher,
    parametrized_packet: RPacket,
    radial_1d_geometry: NumbaRadial1DGeometry,
    classic_opacity_state: OpacityStateNumba,
    bulk_estimators: EstimatorsBulk,
    line_estimators: EstimatorsLine,
    recording_tracker: _RecordingTracker,
    montecarlo_configuration: MonteCarloConfiguration,
) -> None:
    """Claim: The packet loop preserves the emitted line's directional cursors.

    Regime: Numba-disabled propagation in a negative-gradient shell, with a
    mocked line event, emission in ID 0, and subsequent boundary events.
    Verification: The next trace receives the exact cursor pair (1, -1)
    assigned by the emission handler, and the packet escapes. This verifies
    orchestration; physical crossings and emission are tested separately.
    """
    nonhomologous_propagation.packet_propagation(
        parametrized_packet,
        radial_1d_geometry,
        0.0,
        classic_opacity_state,
        bulk_estimators,
        line_estimators,
        recording_tracker,
        montecarlo_configuration,
    )

    assert parametrized_packet.status == PacketStatus.EMITTED
    assert emitted_line_cursors.trace_inputs[1] == (1, -1)


@pytest.mark.parametrize(
    "opacity_state_args",
    [{"line_list_nu": [1.0e15], "tau_sobolev": np.array([[1.0e6]])}],
    indirect=True,
)
def test_nonhomologous_packet_loop_retains_coherent_line_cursor(
    scatter_packet_loop: Callable[..., None],
    line_search_packet: RPacket,
    line_search_geometry: NumbaRadial1DGeometry,
    classic_opacity_state: OpacityStateNumba,
    bulk_estimators: EstimatorsBulk,
    line_estimators: EstimatorsLine,
    recording_tracker: _RecordingTracker,
    montecarlo_configuration: MonteCarloConfiguration,
    set_seed_fixture: Callable[[int], None],
) -> None:
    """The packet loop preserves a crossed line's history in scatter mode.

    Claim: A monotonic outgoing ray cannot replay its emitting line.
    Regime: One nonhomologous shell with positive projected gradient.
    Verification: An analytic first resonance and one line estimator entry.
    """
    packet = line_search_packet
    geometry = line_search_geometry
    initial_radius = packet.r
    initial_mu = packet.mu
    first_resonance_distance = 1.0e12
    resonance_radius = np.sqrt(
        initial_radius**2
        + first_resonance_distance**2
        + 2.0 * initial_radius * first_resonance_distance * initial_mu
    )
    resonance_mu = (
        initial_radius * initial_mu + first_resonance_distance
    ) / resonance_radius
    initial_velocity = geometry.get_velocity(initial_radius, 0)
    resonance_velocity = geometry.get_velocity(resonance_radius, 0)
    line_frequency = classic_opacity_state.line_list_nu[0]
    lab_frequency = line_frequency / (
        1.0 - resonance_velocity * resonance_mu / const.c.cgs.value
    )
    initial_comoving_frequency = lab_frequency * (
        1.0 - initial_velocity * initial_mu / const.c.cgs.value
    )
    packet.nu = initial_comoving_frequency
    classic_opacity_state.electron_density[:] = 1.0e-6
    set_seed_fixture(1963)

    scatter_packet_loop(
        packet,
        geometry,
        0.0,
        classic_opacity_state,
        bulk_estimators,
        line_estimators,
        recording_tracker,
        montecarlo_configuration,
    )

    assert packet.status == PacketStatus.EMITTED
    assert recording_tracker.events.count(InteractionType.LINE) == 1
    npt.assert_allclose(packet.r, geometry.r_outer[0], rtol=1.0e-12)
    expected_line_energy = (
        line_frequency
        / lab_frequency
        / (1.0 - initial_velocity * initial_mu / const.c.cgs.value)
    )
    npt.assert_allclose(
        line_estimators.energy_deposition_line_rate[0, 0],
        expected_line_energy,
        rtol=1.0e-12,
    )


@pytest.mark.xfail(
    reason="Need to update for new mode architecture with correct parameters"
)
# TODO: Update test to provide all required parameters for packet_propagation
def test_verysimple_single_packet_loop(
    verysimple_numba_homologous_radial_1d_geometry,
    verysimple_time_explosion,
    verysimple_opacity_state,
    verysimple_estimators,
    verysimple_vpacket_collection,
    verysimple_packet_collection,
):
    pytest.skip("Test needs to be updated for new mode architecture")
    numba_radial_1d_geometry = verysimple_numba_homologous_radial_1d_geometry
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
    # estimators_bulk, estimators_line, rpacket_tracker,
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
