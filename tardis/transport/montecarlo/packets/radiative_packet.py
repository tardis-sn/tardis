import numba as nb
import numpy as np
from enum import IntEnum
from numba import njit, objmode
from numba.experimental import jitclass

from tardis.transport.frame_transformations import get_doppler_factor
from tardis.transport.montecarlo import njit_dict_no_parallel


class InteractionType(IntEnum):
    """
    Enumeration for different types of packet interactions.

    The values are chosen to enable binary masking operations, allowing
    multiple interaction types to be combined using bitwise operations.
    For example, a packet that experienced both line interaction and
    electron scattering could be represented as LINE | ESCATTERING = 6.

    NO_INTERACTION : int
        Initial state, packet has not yet interacted (value = -1)
    BOUNDARY : int
        Packet crossed a shell boundary (value = 1)
    LINE : int
        Packet underwent line interaction/absorption (value = 2)
    ESCATTERING : int
        Packet underwent electron scattering (value = 4)
    CONTINUUM_PROCESS : int
        Packet underwent continuum process interaction (value = 8)
    """
    NO_INTERACTION = -1
    BOUNDARY = 1
    LINE = 2
    ESCATTERING = 4
    CONTINUUM_PROCESS = 8


class PacketStatus(IntEnum):
    IN_PROCESS = 0
    EMITTED = 1
    REABSORBED = 2
    ADIABATIC_COOLING = 4


@jitclass
class RPacket:
    r: nb.float64  # type: ignore[misc]
    mu: nb.float64  # type: ignore[misc]
    nu: nb.float64  # type: ignore[misc]
    energy: nb.float64  # type: ignore[misc]
    next_line_id: nb.int64  # type: ignore[misc]
    current_shell_id: nb.int64  # type: ignore[misc]
    status: nb.int64  # type: ignore[misc]
    seed: nb.int64  # type: ignore[misc]
    index: nb.int64  # type: ignore[misc]

    def __init__(
        self,
        r: float,
        mu: float,
        nu: float,
        energy: float,
        seed: int,
        index: int = 0,
    ) -> None:
        """
        Initialize radiative packet for Monte Carlo transport.

        Parameters
        ----------
        r : float
            Initial radius [cm].
        mu : float
            Initial directional cosine.
        nu : float
            Initial frequency [Hz].
        energy : float
            Initial energy. Energy units are scaled with time_simulation.
            Adds all up to 1 in a single run.
        seed : int
            Random number seed.
        index : int, optional
            Packet index, by default 0.
        """
        self.r = r
        self.mu = mu
        self.nu = nu
        self.energy = energy
        self.current_shell_id = 0
        self.status = PacketStatus.IN_PROCESS
        self.seed = seed
        self.index = index

    def initialize_line_id(
        self, opacity_state, time_explosion, enable_full_relativity
    ):
        inverse_line_list_nu = opacity_state.line_list_nu[::-1]
        doppler_factor = get_doppler_factor(
            self.r, self.mu, time_explosion, enable_full_relativity
        )
        comov_nu = self.nu * doppler_factor
        next_line_id = len(opacity_state.line_list_nu) - np.searchsorted(
            inverse_line_list_nu, comov_nu
        )
        if next_line_id == len(opacity_state.line_list_nu):
            next_line_id -= 1
        self.next_line_id = next_line_id


@njit(**njit_dict_no_parallel)
def print_r_packet_properties(r_packet):
    """
    Print all packet information.

    Parameters
    ----------
    r_packet : RPacket
        RPacket object.
    """
    print("-" * 80)
    print("R-Packet information:")
    with objmode:
        print("r =", str(r_packet.r))
        print("mu =", str(r_packet.mu))
        print("nu =", str(r_packet.nu))
        print("energy =", str(r_packet.energy))
        print("next_line_id =", str(r_packet.next_line_id))
        print("current_shell_id =", str(r_packet.current_shell_id))
        print("status =", str(r_packet.status))
        print("seed =", str(r_packet.seed))
        print("index =", str(r_packet.index))
        print("last_interaction_type =", str(r_packet.last_interaction_type))
        print("last_interaction_in_nu =", str(r_packet.last_interaction_in_nu))
        print("last_interaction_in_r =", str(r_packet.last_interaction_in_r))
        print("last_line_interaction_in_id =", str(r_packet.last_line_interaction_in_id))
        print("last_line_interaction_out_id =", str(r_packet.last_line_interaction_out_id))
        print("last_line_interaction_shell_id =", str(r_packet.last_line_interaction_shell_id))
    print("-" * 80)
