import numpy as np
import numpy.testing as npt
import pytest

import tardis.transport.montecarlo.packets.virtual_packet as virtual_packet
from tardis import constants as const
from tardis.transport.frame_transformations import (
    get_doppler_factor,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    VPacketCollection,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
)
from tardis.transport.montecarlo.packets.trackers.tracker_last_interaction import (
    TrackerLastInteraction,
)

C_SPEED_OF_LIGHT = const.c.to("cm/s").value


@pytest.fixture(scope="function")
def v_packet():
    return virtual_packet.VPacket(
        r=7.5e14,
        nu=4e15,
        mu=0.3,
        energy=0.9,
        current_shell_id=0,
        next_line_id=0,
        index=0,
    )


def v_packet_initialize_line_id(v_packet, opacity_state, time_explosion):
    inverse_line_list_nu = opacity_state.line_list_nu[::-1]
    velocity = v_packet.r / time_explosion
    doppler_factor = get_doppler_factor(velocity, v_packet.mu, False)
    comov_nu = v_packet.nu * doppler_factor
    next_line_id = len(opacity_state.line_list_nu) - np.searchsorted(
        inverse_line_list_nu, comov_nu
    )
    v_packet.next_line_id = next_line_id


def test_store_vpacket_preserves_last_line_interaction_metadata():
    vpacket_collection = VPacketCollection(
        source_rpacket_index=0,
        spectrum_frequency_grid=np.array([1.0, 2.0]),
        number_of_vpackets=1,
        v_packet_spawn_start_frequency=0.0,
        v_packet_spawn_end_frequency=np.inf,
        temporary_v_packet_bins=1,
    )
    v_packet = virtual_packet.VPacket(
        r=5.0,
        nu=2.0,
        mu=0.3,
        energy=1.0,
        current_shell_id=4,
        next_line_id=3,
    )
    rpacket_tracker = TrackerLastInteraction()
    rpacket_tracker.before_nu = 7.0
    rpacket_tracker.radius = 8.0
    rpacket_tracker.interaction_type = InteractionType.LINE
    rpacket_tracker.interaction_line_absorb_id = 11
    rpacket_tracker.interaction_line_emit_id = 12
    rpacket_tracker.shell_id = 4

    vpacket_collection.add_packet_from_tracker(
        v_packet.nu,
        v_packet.energy,
        v_packet.mu,
        v_packet.r,
        rpacket_tracker,
    )

    npt.assert_allclose(vpacket_collection.last_interaction_in_nu, [7.0])
    npt.assert_allclose(vpacket_collection.last_interaction_in_r, [8.0])
    npt.assert_array_equal(
        vpacket_collection.last_interaction_type, [InteractionType.LINE]
    )
    npt.assert_array_equal(vpacket_collection.last_interaction_in_id, [11])
    npt.assert_array_equal(vpacket_collection.last_interaction_out_id, [12])
    npt.assert_array_equal(vpacket_collection.last_interaction_shell_id, [4])


def test_store_vpacket_preserves_no_interaction_defaults():
    vpacket_collection = VPacketCollection(
        source_rpacket_index=0,
        spectrum_frequency_grid=np.array([1.0, 2.0]),
        number_of_vpackets=1,
        v_packet_spawn_start_frequency=0.0,
        v_packet_spawn_end_frequency=np.inf,
        temporary_v_packet_bins=1,
    )
    rpacket_tracker = TrackerLastInteraction()

    vpacket_collection.add_packet_from_tracker(
        2.0,
        1.0,
        0.3,
        5.0,
        rpacket_tracker,
    )

    assert np.isnan(vpacket_collection.last_interaction_in_nu[0])
    assert np.isnan(vpacket_collection.last_interaction_in_r[0])
    npt.assert_array_equal(
        vpacket_collection.last_interaction_type,
        [InteractionType.NO_INTERACTION],
    )
    npt.assert_array_equal(vpacket_collection.last_interaction_in_id, [-1])
    npt.assert_array_equal(vpacket_collection.last_interaction_out_id, [-1])
    npt.assert_array_equal(vpacket_collection.last_interaction_shell_id, [-1])


def test_store_vpacket_preserves_escattering_metadata():
    vpacket_collection = VPacketCollection(
        source_rpacket_index=0,
        spectrum_frequency_grid=np.array([1.0, 2.0]),
        number_of_vpackets=1,
        v_packet_spawn_start_frequency=0.0,
        v_packet_spawn_end_frequency=np.inf,
        temporary_v_packet_bins=1,
    )
    rpacket_tracker = TrackerLastInteraction()
    rpacket_tracker.before_nu = 7.0
    rpacket_tracker.radius = 8.0
    rpacket_tracker.interaction_type = InteractionType.ESCATTERING
    rpacket_tracker.shell_id = 4

    vpacket_collection.add_packet_from_tracker(
        2.0,
        1.0,
        0.3,
        5.0,
        rpacket_tracker,
    )

    npt.assert_allclose(vpacket_collection.last_interaction_in_nu, [7.0])
    npt.assert_allclose(vpacket_collection.last_interaction_in_r, [8.0])
    npt.assert_array_equal(
        vpacket_collection.last_interaction_type,
        [InteractionType.ESCATTERING],
    )
    npt.assert_array_equal(vpacket_collection.last_interaction_in_id, [-1])
    npt.assert_array_equal(vpacket_collection.last_interaction_out_id, [-1])
    npt.assert_array_equal(vpacket_collection.last_interaction_shell_id, [4])


def test_trace_vpacket_within_shell(
    v_packet,
    verysimple_numba_homologous_radial_1d_geometry,
    verysimple_time_explosion,
    verysimple_opacity_state,
):
    # Give the vpacket a reasonable line ID
    v_packet_initialize_line_id(
        v_packet, verysimple_opacity_state, verysimple_time_explosion
    )

    (
        tau_trace_combined,
        distance_boundary,
        delta_shell,
    ) = virtual_packet.trace_vpacket_within_shell(
        v_packet,
        verysimple_numba_homologous_radial_1d_geometry,
        verysimple_time_explosion,
        verysimple_opacity_state,
        enable_full_relativity=False,
    )

    npt.assert_almost_equal(tau_trace_combined, 8164850.891288479)
    # changed from almost equal to allclose. Now seems to work.
    npt.assert_allclose(distance_boundary, 843684056256104.1)
    assert delta_shell == 1


def test_trace_vpacket(
    v_packet,
    verysimple_numba_homologous_radial_1d_geometry,
    verysimple_time_explosion,
    verysimple_opacity_state,
):
    # Set seed because of RNG in trace_vpacket
    np.random.seed(1)

    # Give the vpacket a reasonable line ID
    v_packet_initialize_line_id(
        v_packet, verysimple_opacity_state, verysimple_time_explosion
    )

    tau_trace_combined = virtual_packet.trace_vpacket(
        v_packet,
        verysimple_numba_homologous_radial_1d_geometry,
        verysimple_time_explosion,
        verysimple_opacity_state,
        10.0,
        0.0,
        enable_full_relativity=False,
    )

    npt.assert_almost_equal(tau_trace_combined, 8164850.891288479)
    # change from almost_equal to allclose. Now seems to work.
    npt.assert_allclose(v_packet.r, 1286064000000000.0)
    npt.assert_almost_equal(v_packet.nu, 4.0e15)
    npt.assert_almost_equal(v_packet.energy, 0.0)
    npt.assert_almost_equal(v_packet.mu, 0.8309726858508629)
    assert v_packet.next_line_id == 2773
    assert v_packet.current_shell_id == 1


# NEEDS TO TEST VPACKET COLLECTION OVERFLOW
@pytest.mark.xfail(reason="Needs to be implemented")
def test_trace_vpacket_volley(
    packet,
    verysimple_packet_collection,
    verysimple_3vpacket_collection,
    verysimple_numba_homologous_radial_1d_geometry,
    verysimple_time_explosion,
    verysimple_opacity_state,
):
    # Set seed because of RNG in trace_vpacket
    np.random.seed(1)

    packet.initialize_line_id(
        verysimple_opacity_state, verysimple_time_explosion
    )
    rpacket_tracker = TrackerLastInteraction()

    virtual_packet.trace_vpacket_volley(
        packet,
        rpacket_tracker,
        verysimple_3vpacket_collection,
        verysimple_numba_homologous_radial_1d_geometry,
        verysimple_time_explosion,
        verysimple_opacity_state,
        enable_full_relativity=False,
        tau_russian=10.0,
        survival_probability=0.0,
    )


@pytest.fixture(scope="function")
def broken_packet():
    return virtual_packet.VPacket(
        r=1286064000000000.0,
        nu=1660428912896553.2,
        mu=0.4916053094346575,
        energy=2.474533071386993e-07,
        index=3,
        current_shell_id=0,
        next_line_id=5495,
    )


def test_trace_bad_vpacket(
    broken_packet,
    verysimple_numba_homologous_radial_1d_geometry,
    verysimple_time_explosion,
    verysimple_opacity_state,
):
    virtual_packet.trace_vpacket(
        broken_packet,
        verysimple_numba_homologous_radial_1d_geometry,
        verysimple_time_explosion,
        verysimple_opacity_state,
        10.0,
        0.0,
        enable_full_relativity=False,
    )
