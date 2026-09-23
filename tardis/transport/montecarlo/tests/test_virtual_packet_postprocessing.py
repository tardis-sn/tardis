import numpy as np
import numpy.testing as npt
import pytest

import tardis.transport.montecarlo.configuration.montecarlo_globals as montecarlo_globals
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.modes.montecarlo_transport import (
    calculate_virtual_packet_spectrum,
)
from tardis.transport.montecarlo.modes.nonhomologous.virtual_packet import (
    trace_vpacket_volley,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    PacketCollection,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    RPacket,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full_util import (
    generate_tracker_full_list,
)


def test_vpacket_postprocessing_uses_postinteraction_tracker_state(
    opacity_state_args: tuple,
    montecarlo_configuration: MonteCarloConfiguration,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify tracked post-interaction states source the virtual spectrum.

    Claim: A tracked line interaction contributes its post-interaction packet
    state to the virtual spectrum.
    Regime: Static transparent ejecta, with only the post-interaction frequency
    inside the virtual-packet spawn range.
    Verification: The virtual energy equals the analytic off-photosphere
    angular weight summed over the volley.
    """
    monkeypatch.setattr(
        montecarlo_globals, "CONTINUUM_PROCESSES_ENABLED", False
    )
    opacity_args = list(opacity_state_args)
    opacity_args[0] = np.zeros(2)
    opacity_args[3] = np.zeros((2, 2))
    opacity_state = OpacityStateNumba(*opacity_args)
    geometry = NumbaRadial1DGeometry(
        np.array([1.0, 2.0]),
        np.array([2.0, 10.0]),
        np.zeros(2),
        np.zeros(2),
    )
    spectrum_frequency_grid = np.array([1.0e14, 2.0e14, 3.0e14])
    packet_collection = PacketCollection(
        np.array([1.0]),
        np.array([1.0e14]),
        np.array([0.5]),
        np.array([1.0]),
        np.array([1963], dtype=np.int64),
        1.0,
    )
    trackers = generate_tracker_full_list(1, 3)
    tracker = trackers[0]
    r_packet = RPacket(1.0, 0.5, 1.0e14, 1.0, 1963)
    r_packet.next_line_id = 0
    tracker.track_boundary_event(r_packet, -1, 0)

    r_packet.r = 4.0
    r_packet.current_shell_id = 1
    r_packet.energy = 0.25
    tracker.track_line_interaction_before(r_packet)
    r_packet.nu = 2.0e14
    r_packet.mu = 0.2
    r_packet.energy = 2.0
    r_packet.next_line_id = 1
    tracker.track_line_interaction_after(r_packet)
    tracker.track_boundary_event(r_packet, 1, 2)
    tracker.finalize()

    montecarlo_configuration.ENABLE_FULL_RELATIVITY = False
    montecarlo_configuration.ENABLE_VPACKET_TRACKING = True
    montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY = 1.5e14
    montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY = 2.5e14
    montecarlo_configuration.TEMPORARY_V_PACKET_BINS = 4

    energy_histogram, vpacket_tracker = calculate_virtual_packet_spectrum(
        packet_collection,
        geometry,
        0.0,
        opacity_state,
        montecarlo_configuration,
        spectrum_frequency_grid,
        trackers,
        4,
        trace_vpacket_volley,
    )

    mu_min = -np.sqrt(1.0 - (geometry.r_inner[0] / r_packet.r) ** 2)
    expected_energy = r_packet.energy * (1.0 - mu_min) / 2.0
    npt.assert_allclose(
        energy_histogram,
        np.array([0.0, expected_energy, 0.0]),
        rtol=1.0e-14,
        atol=0.0,
    )
    npt.assert_array_equal(
        vpacket_tracker.last_interaction_type,
        np.full(4, InteractionType.LINE),
    )
    npt.assert_array_equal(
        vpacket_tracker.last_interaction_in_nu,
        np.full(4, 1.0e14),
    )
