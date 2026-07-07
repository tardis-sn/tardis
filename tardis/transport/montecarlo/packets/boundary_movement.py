"""Shared shell-boundary movement for radiative packets."""

from numba import njit

from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.packets.radiative_packet import (
    PacketStatus,
    RPacket,
)


@njit(**njit_dict_no_parallel)
def move_packet_across_shell_boundary(
    packet: RPacket, delta_shell: int, no_of_shells: int
) -> None:
    """
    Move a packet across a shell boundary and update its transport status.

    Parameters
    ----------
    packet : RPacket
        Radiative packet object.
    delta_shell : int
        Change in shell index, positive outward and negative inward.
    no_of_shells : int
        Number of shells in the transport grid.
    """
    next_shell_id = packet.current_shell_id + delta_shell

    if next_shell_id >= no_of_shells:
        packet.status = PacketStatus.EMITTED
    elif next_shell_id < 0:
        packet.status = PacketStatus.REABSORBED
    else:
        packet.current_shell_id = next_shell_id

