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
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
    init_estimators_bulk,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
    init_estimators_line,
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

RTOL = 1.0e-12


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
        self.trace_inputs: list[tuple[object, object]] = []
        self.absorbed_line_id: int | None = None

    def trace_packet(
        self, *args: object, **kwargs: object
    ) -> tuple[object, ...]:
        assert self.interaction_queue
        self.trace_inputs.append((args[-2], args[-1]))
        interaction = self.interaction_queue.pop(0)
        if (
            self.absorbed_line_id is not None
            and interaction[1] == InteractionType.LINE
        ):
            args[0].next_line_id = self.absorbed_line_id
        return interaction

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
    # Force the nonhomologous loop through one selected interaction, then an
    # outer boundary, so this characterizes dispatch bookkeeping directly.
    interactions = [(1.0e12, first_interaction, 1, 0, 0)]
    if first_interaction != InteractionType.BOUNDARY:
        interactions.append((1.0e12, InteractionType.BOUNDARY, 1, 0, 0))
    interactions.append((1.0e12, InteractionType.BOUNDARY, 1, 0, 0))
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
def test_nonhomologous_packet_loop_reseeds_after_different_line_emission(
    python_numba_disabled: None,
    patch_common_classic_hooks: _CommonPacketPropagationPatcher,
    parametrized_packet: RPacket,
    radial_1d_geometry: NumbaRadial1DGeometry,
    classic_opacity_state: OpacityStateNumba,
    bulk_estimators: EstimatorsBulk,
    line_estimators: EstimatorsLine,
    recording_tracker: _RecordingTracker,
    montecarlo_configuration: MonteCarloConfiguration,
) -> None:
    """A different emitted line seeds a fresh directional cursor.

    This test checks packet-loop state transfer; direct transport tests check
    physical line positions and estimator contributions.
    """
    interactions = [
        (1.0e12, InteractionType.LINE, 1, 2, 1),
        (1.0e12, InteractionType.BOUNDARY, 1, 1, 1),
        (1.0e12, InteractionType.BOUNDARY, 1, 1, 1),
    ]
    patch_common_classic_hooks(nonhomologous_propagation, interactions)
    patch_common_classic_hooks.absorbed_line_id = 1

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
    assert patch_common_classic_hooks.trace_inputs[1] == (1, 1)


@pytest.mark.parametrize(
    "opacity_state_args",
    [{"line_list_nu": [1.0e15], "tau_sobolev": np.array([[1.0e6]])}],
    indirect=True,
)
def test_nonhomologous_packet_loop_retains_coherent_line_cursor(
    monkeypatch: pytest.MonkeyPatch,
    classic_opacity_state: OpacityStateNumba,
    recording_tracker: _RecordingTracker,
    montecarlo_configuration: MonteCarloConfiguration,
    set_seed_fixture: Callable[[int], None],
) -> None:
    """The packet loop carries a crossed line past coherent emission.

    Claim: A monotonic outgoing ray cannot replay its emitting line.
    Regime: One nonhomologous shell with positive projected gradient.
    Verification: An analytic first resonance and one line estimator entry.
    """
    geometry = NumbaRadial1DGeometry(
        np.array([1.0e14]),
        np.array([2.0e14]),
        np.array([1.0e8]),
        np.array([1.8e8]),
    )
    initial_radius = 1.25e14
    initial_mu = 0.1
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
    packet = RPacket(
        r=initial_radius,
        mu=initial_mu,
        nu=initial_comoving_frequency,
        energy=1.0,
        seed=1963,
    )
    classic_opacity_state.electron_density[:] = 1.0e-6
    bulk_estimators = init_estimators_bulk(1)
    line_estimators = init_estimators_line((1, 1))
    recording_tracker.max_line_events = 1
    montecarlo_configuration.LINE_INTERACTION_TYPE = 0
    montecarlo_configuration.DISABLE_LINE_SCATTERING = False
    set_seed_fixture(1963)

    # Keep the real event physics and tracer; control only the angle draw.
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
    packet_loop = getattr(
        nonhomologous_propagation.packet_propagation,
        "py_func",
        nonhomologous_propagation.packet_propagation,
    )
    packet_loop(
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
