from enum import IntEnum

import numpy as np
import pandas as pd
from numba import int64, float64, njit, objmode
from numba.experimental import jitclass

from tardis.montecarlo.montecarlo_numba import (
    njit_dict_no_parallel,
)
from tardis.transport.frame_transformations import (
    get_doppler_factor,
)
from tardis.montecarlo.montecarlo_numba import numba_config as nc
from tardis.montecarlo.montecarlo_numba import njit_dict_no_parallel


class InteractionType(IntEnum):
    BOUNDARY = 1
    LINE = 2
    ESCATTERING = 4
    CONTINUUM_PROCESS = 8


class PacketStatus(IntEnum):
    IN_PROCESS = 0
    EMITTED = 1
    REABSORBED = 2
    ADIABATIC_COOLING = 4


rpacket_spec = [
    ("r", float64),
    ("mu", float64),
    ("nu", float64),
    ("energy", float64),
    ("next_line_id", int64),
    ("current_shell_id", int64),
    ("status", int64),
    ("seed", int64),
    ("index", int64),
    ("last_interaction_type", int64),
    ("last_interaction_in_nu", float64),
    ("last_line_interaction_in_id", int64),
    ("last_line_interaction_out_id", int64),
]


@jitclass(rpacket_spec)
class RPacket(object):
    def __init__(self, r, mu, nu, energy, seed, index=0):
        self.r = r
        self.mu = mu
        self.nu = nu
        self.energy = energy
        self.current_shell_id = 0
        self.status = PacketStatus.IN_PROCESS
        self.seed = seed
        self.index = index
        self.last_interaction_type = -1
        self.last_interaction_in_nu = 0.0
        self.last_line_interaction_in_id = -1
        self.last_line_interaction_out_id = -1

    def initialize_line_id(self, numba_plasma, numba_model):
        inverse_line_list_nu = numba_plasma.line_list_nu[::-1]
        doppler_factor = get_doppler_factor(
            self.r, self.mu, numba_model.time_explosion
        )
        comov_nu = self.nu * doppler_factor
        next_line_id = len(numba_plasma.line_list_nu) - np.searchsorted(
            inverse_line_list_nu, comov_nu
        )
        if next_line_id == len(numba_plasma.line_list_nu):
            next_line_id -= 1
        self.next_line_id = next_line_id


@njit(**njit_dict_no_parallel)
def print_r_packet_properties(r_packet):
    """
    Print all packet information

    Parameters
    ----------
    r_packet : RPacket
        RPacket object
    """
    print("-" * 80)
    print("R-Packet information:")
    with objmode:
        for r_packet_attribute_name, _ in rpacket_spec:
            print(
                r_packet_attribute_name,
                "=",
                str(getattr(r_packet, r_packet_attribute_name)),
            )
    print("-" * 80)


def rpacket_trackers_to_dataframe(rpacket_trackers):
    """Generates a dataframe from the rpacket_trackers list of RPacketCollection Objects.

    Parameters
    ----------
    rpacket_trackers : numba.typed.typedlist.List
        list of individual RPacketCollection class objects

    Returns
    -------
    pandas.core.frame.DataFrame
        Dataframe containing properties of RPackets
    """
    rtracker_dict_list = []
    index_array = []
    step_array = []
    for rpacket in rpacket_trackers:
        for rpacket_step_no in range(len(rpacket.r)):
            index_array.append(rpacket.index)
            step_array.append(rpacket_step_no)
            rtracker_dict_list.append(
                [
                    rpacket.status[rpacket_step_no],
                    rpacket.seed,
                    rpacket.r[rpacket_step_no],
                    rpacket.nu[rpacket_step_no],
                    rpacket.mu[rpacket_step_no],
                    rpacket.energy[rpacket_step_no],
                    rpacket.shell_id[rpacket_step_no],
                ]
            )
    index_array = [np.array(index_array), np.array(step_array)]
    rpacket_df = pd.DataFrame(
        rtracker_dict_list,
        index=index_array,
        columns=["status", "seed", "r", "nu", "mu", "energy", "shell_id"],
    )
    rpacket_df.index.names = ["index", "step"]
    return rpacket_df
