import pytest

import numpy.testing as npt

from tardis.montecarlo.montecarlo_numba import RPacket
from tardis.montecarlo.montecarlo_numba.single_packet_loop import (
    single_packet_loop,
)


@pytest.mark.xfail(reason="Need to fix estimator differences across runs")
# TODO set RNG consistently
def test_verysimple_single_packet_loop(
    verysimple_numba_radial_1d_geometry,
    verysimple_numba_model,
    verysimple_numba_plasma,
    verysimple_estimators,
    verysimple_vpacket_collection,
    verysimple_packet_collection,
):
    numba_radial_1d_geometry = verysimple_numba_radial_1d_geometry
    packet_collection = verysimple_packet_collection
    vpacket_collection = verysimple_vpacket_collection
    numba_model = verysimple_numba_model
    numba_plasma = verysimple_numba_plasma
    numba_estimators = verysimple_estimators

    i = 0
    r_packet = RPacket(
        numba_radial_1d_geometry.r_inner[0],
        packet_collection.packets_input_mu[i],
        packet_collection.packets_input_nu[i],
        packet_collection.packets_input_energy[i],
        i,
    )
    single_packet_loop(
        r_packet,
        numba_radial_1d_geometry,
        numba_model,
        numba_plasma,
        numba_estimators,
        vpacket_collection,
    )

    npt.assert_almost_equal(r_packet.nu, 1053057938883272.8)
    npt.assert_almost_equal(r_packet.mu, 0.9611146425440562)
    npt.assert_almost_equal(r_packet.energy, 0.10327717505563379)


@pytest.mark.xfail(reason="To be implemented")
def test_set_packet_props_partial_relativity():
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_set_packet_props_full_relativity():
    assert False
