import numpy as np
from numba import njit

from tardis.montecarlo.montecarlo_numba import (
    njit_dict_no_parallel,
)

from tardis.montecarlo.montecarlo_numba.numba_config import ENABLE_NONHOMOLOGOUS_EXPANSION

@njit(**njit_dict_no_parallel)
def velocity(r_packet, shell_id, numba_model):
    """
    Velocity at radius r

    Parameters
    ----------
    r_packet: RPacket
    shell_id: int
    numba_model: NumbaModel

    Returns
    -----------
    v: float, current velocity
    """
    v_inner = numba_model.v_inner[shell_id]
    v_outer = numba_model.v_outer[shell_id]
    r_inner = numba_model.r_inner[shell_id]
    r_outer = numba_model.r_outer[shell_id]
    r = r_packet.r
    frac = (v_outer-v_inner)/(r_outer-r_inner)
    return v_inner + frac * (r - r_inner)

@njit(**njit_dict_no_parallel)
def dvdr(shell_id, numba_model):
    """
    dv/dr for the current shell

    Parameters
    ----------
    shell_id: int
    numba_model: NumbaModel

    Returns
    -----------
    dvdr: float, dv/dr of the current shell
    """
    v_inner = numba_model.v_inner[shell_id]
    v_outer = numba_model.v_outer[shell_id]
    r_inner = numba_model.r_inner[shell_id]
    r_outer = numba_model.r_outer[shell_id]
    return (v_outer-v_inner)/(r_outer-r_inner)


def tau_sobolev_factor(r_packet, shell_id, numba_model, numba_plasma):
    if ENABLE_NONHOMOLOGOUS_EXPANSION:
        v = velocity(r_packet, shell_id, numba_model)
        r = r_packet.r
        dvdr = dvdr(shell_id, numba_model)
        mu = r_packet.mu
        factor = 1.0/((1 - mu * mu)*v/r + mu * mu * dvdr)


    else:
        factor = np.ones(numba_plasma.tau_sobolev.shape)
    return factor * numba_plasma.tau_sobolev

