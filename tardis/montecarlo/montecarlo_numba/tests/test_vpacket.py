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

from numpy.testing import (
        assert_equal,
        assert_almost_equal,
        assert_array_equal,
        assert_allclose
        )


@pytest.fixture(scope="function")
def v_packet():
    return vpacket.VPacket(
        r = 7.5e14,
        nu = 0.4,
        mu = 0.3,
        energy = 0.9,
        current_shell_id = 0,
        next_line_id = 0,
        index = 0,
        is_close_line = 0
    )

@pytest.fixture(scope="function")
def model():
    return numba_interface.NumbaModel(
        r_inner = np.array([6.912e14, 8.64e14], dtype=np.float64),
        r_outer = np.array([8.64e14, 1.0368e15], dtype=np.float64),
        time_explosion = 5.2e7
    )

@pytest.fixture(scope="function")
def plasma():
    return numba_interface.NumbaPlasma(
        electron_density=1.0e9*np.ones(2, dtype=np.float64),
        line_list_nu=np.array([
            1.26318289e+16,
            1.26318289e+16,
            1.23357675e+16,
            1.23357675e+16,
            1.16961598e+16], dtype=np.float64),
        tau_sobolev=np.ones((2, 1000), dtype=np.float64),
        transition_probabilities=np.zeros((2, 2), dtype=np.float64),
        line2macro_level_upper=np.zeros(2, dtype=np.int64),
        macro_block_references=np.zeros(2, dtype=np.int64),
        transition_type=np.zeros(2, dtype=np.int64),
        destination_level_id=np.zeros(2, dtype=np.int64),
        transition_line_id=np.zeros(2, dtype=np.int64)
    )


@pytest.mark.xfail(reason='To be implemented')
def test_trace_vpacket_within_shell(v_packet, model, plasma):
    #tau_trace_combined, distance_boundary, delta_shell = vpacket.trace_vpacket_within_shell(
    #                                                                                v_packet,
    #                                                                                model,
    #                                                                                plasma)
    assert False


@pytest.mark.xfail(reason='To be implemented')
def test_trace_vpacket():
    assert False

@pytest.mark.xfail(reason='To be implemented')
def test_trace_vpacket_volley():
    assert False

