import numpy.testing as npt
import pytest

from tardis.transport.montecarlo import RPacket
from tardis.transport.montecarlo.modes.classic.packet_propagation import (
    packet_propagation,
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
