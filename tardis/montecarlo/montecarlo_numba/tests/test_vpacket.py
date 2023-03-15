import os
import pytest
import numpy as np
import pandas as pd
import tardis.montecarlo.montecarlo_numba.r_packet as r_packet
import tardis.montecarlo.montecarlo_numba.vpacket as vpacket
import tardis.montecarlo.montecarlo_configuration as mc
import tardis.montecarlo.montecarlo_numba.numba_interface as numba_interface
from tardis import constants as const
from tardis.montecarlo.montecarlo_numba import macro_atom

from tardis.transport.frame_transformations import (
    get_doppler_factor,
)

C_SPEED_OF_LIGHT = const.c.to("cm/s").value

import numpy.testing as npt


@pytest.fixture(scope="function")
def v_packet():
    return vpacket.VPacket(
        r=7.5e14,
        nu=4e15,
        mu=0.3,
        energy=0.9,
        current_shell_id=0,
        next_line_id=0,
        index=0,
    )


def v_packet_initialize_line_id(v_packet, numba_plasma, numba_model):
    inverse_line_list_nu = numba_plasma.line_list_nu[::-1]
    doppler_factor = get_doppler_factor(
        v_packet.r, v_packet.mu, numba_model.time_explosion
    )
    comov_nu = v_packet.nu * doppler_factor
    next_line_id = len(numba_plasma.line_list_nu) - np.searchsorted(
        inverse_line_list_nu, comov_nu
    )
    v_packet.next_line_id = next_line_id


def test_trace_vpacket_within_shell(
    v_packet,
    verysimple_numba_radial_1d_geometry,
    verysimple_numba_model,
    verysimple_numba_plasma,
):
    # Give the vpacket a reasonable line ID
    v_packet_initialize_line_id(
        v_packet, verysimple_numba_plasma, verysimple_numba_model
    )

    (
        tau_trace_combined,
        distance_boundary,
        delta_shell,
    ) = vpacket.trace_vpacket_within_shell(
        v_packet,
        verysimple_numba_radial_1d_geometry,
        verysimple_numba_model,
        verysimple_numba_plasma,
    )

    npt.assert_almost_equal(tau_trace_combined, 8164850.891288479)
    npt.assert_almost_equal(distance_boundary, 843684056256104.1)
    assert delta_shell == 1


def test_trace_vpacket(
    v_packet,
    verysimple_numba_radial_1d_geometry,
    verysimple_numba_model,
    verysimple_numba_plasma,
):
    # Set seed because of RNG in trace_vpacket
    np.random.seed(1)

    # Give the vpacket a reasonable line ID
    v_packet_initialize_line_id(
        v_packet, verysimple_numba_plasma, verysimple_numba_model
    )

    tau_trace_combined = vpacket.trace_vpacket(
        v_packet,
        verysimple_numba_radial_1d_geometry,
        verysimple_numba_model,
        verysimple_numba_plasma,
    )

    npt.assert_almost_equal(tau_trace_combined, 8164850.891288479)
    npt.assert_almost_equal(v_packet.r, 1286064000000000.0)
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
    verysimple_numba_radial_1d_geometry,
    verysimple_numba_model,
    verysimple_numba_plasma,
):
    # Set seed because of RNG in trace_vpacket
    np.random.seed(1)

    packet.initialize_line_id(verysimple_numba_plasma, verysimple_numba_model)

    vpacket.trace_vpacket_volley(
        packet,
        verysimple_3vpacket_collection,
        verysimple_numba_radial_1d_geometry,
        verysimple_numba_model,
        verysimple_numba_plasma,
    )


@pytest.fixture(scope="function")
def broken_packet():
    return vpacket.VPacket(
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
    verysimple_numba_radial_1d_geometry,
    verysimple_numba_model,
    verysimple_numba_plasma,
):
    vpacket.trace_vpacket(
        broken_packet,
        verysimple_numba_radial_1d_geometry,
        verysimple_numba_model,
        verysimple_numba_plasma,
    )
