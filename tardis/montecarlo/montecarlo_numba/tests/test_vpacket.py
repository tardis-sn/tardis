import os
import pytest
import numpy as np
import pandas as pd
import tardis.montecarlo.formal_integral as formal_integral
import tardis.montecarlo.montecarlo_numba.r_packet as r_packet
import tardis.montecarlo.montecarlo_numba.vpacket as vpacket
import tardis.montecarlo.montecarlo_configuration as mc
import tardis.montecarlo.montecarlo_numba.numba_interface as numba_interface
from tardis import constants as const
from tardis.montecarlo.montecarlo_numba.numba_interface import Estimators
from tardis.montecarlo.montecarlo_numba import macro_atom
C_SPEED_OF_LIGHT = const.c.to('cm/s').value

import numpy.testing as npt


@pytest.fixture(scope="function")
def v_packet():
    return vpacket.VPacket(
        r = 7.5e14,
        nu = 4e15,
        mu = 0.3,
        energy = 0.9,
        current_shell_id = 0,
        next_line_id = 0,
        index = 0,
        is_close_line = 0
    )

def v_packet_initialize_line_id(v_packet, numba_plasma, numba_model):
        inverse_line_list_nu = numba_plasma.line_list_nu[::-1]
        doppler_factor = r_packet.get_doppler_factor(v_packet.r, v_packet.mu,
                                            numba_model.time_explosion)
        comov_nu = v_packet.nu * doppler_factor
        next_line_id = (len(numba_plasma.line_list_nu) -
                        np.searchsorted(inverse_line_list_nu, comov_nu))
        v_packet.next_line_id = next_line_id

#@pytest.mark.xfail(reason='To be implemented')
def test_trace_vpacket_within_shell(v_packet, verysimple_numba_model, verysimple_numba_plasma):
    #Give the vpacket a reasonable line ID
    v_packet_initialize_line_id(v_packet, verysimple_numba_plasma, verysimple_numba_model)

    tau_trace_combined, distance_boundary, delta_shell = vpacket.trace_vpacket_within_shell(
                                                                                    v_packet,
                                                                                    verysimple_numba_model,
                                                                                    verysimple_numba_plasma)

    npt.assert_almost_equal(tau_trace_combined, 8164850.891288479)
    npt.assert_almost_equal(distance_boundary, 843684056256104.1)
    assert delta_shell == 1


@pytest.mark.xfail(reason='To be implemented')
def test_trace_vpacket():
    assert False

@pytest.mark.xfail(reason='To be implemented')
def test_trace_vpacket_volley():
    assert False

