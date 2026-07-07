"""Packet-level frame transformations for radiative transport."""

from numba import njit

from tardis.transport.frame_transformations import (
    angle_aberration_CMF_to_LF_from_velocity,
    angle_aberration_LF_to_CMF_from_velocity,
)
from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.packets.radiative_packet import RPacket


@njit(**njit_dict_no_parallel)
def transform_packet_lab_to_comoving_frame(
    r_packet: RPacket, geometry, enable_full_relativity: bool
) -> tuple[float, float, float]:
    """
    Transform packet frequency, energy, and angle from lab to comoving frame.

    Parameters
    ----------
    r_packet : RPacket
        Radiative packet object.
    geometry
        Geometry object exposing ``get_velocity(r, shell_id)``.
    enable_full_relativity : bool
        Flag to enable full relativistic calculations.

    Returns
    -------
    tuple[float, float, float]
        Comoving-frame frequency, energy, and angle cosine.
    """
    doppler_factor = geometry.get_doppler_factor(
        r_packet.r,
        r_packet.mu,
        r_packet.current_shell_id,
        enable_full_relativity,
    )
    comoving_mu = r_packet.mu
    if enable_full_relativity:
        velocity = geometry.get_velocity(
            r_packet.r,
            r_packet.current_shell_id,
        )
        comoving_mu = angle_aberration_LF_to_CMF_from_velocity(
            velocity,
            r_packet.mu,
        )
    return (
        r_packet.nu * doppler_factor,
        r_packet.energy * doppler_factor,
        comoving_mu,
    )


@njit(**njit_dict_no_parallel)
def transform_packet_comoving_to_lab_frame(
    r_packet: RPacket,
    geometry,
    comoving_nu: float,
    comoving_energy: float,
    comoving_mu: float,
    enable_full_relativity: bool,
) -> None:
    """
    Transform packet frequency, energy, and angle from comoving to lab frame.

    Parameters
    ----------
    r_packet : RPacket
        Radiative packet object to update in place.
    geometry
        Geometry object exposing ``get_velocity(r, shell_id)``.
    comoving_nu : float
        Comoving-frame packet frequency.
    comoving_energy : float
        Comoving-frame packet energy.
    comoving_mu : float
        Comoving-frame packet propagation angle cosine.
    enable_full_relativity : bool
        Flag to enable full relativistic calculations.
    """
    inverse_doppler_factor = geometry.get_inverse_doppler_factor(
        r_packet.r,
        comoving_mu,
        r_packet.current_shell_id,
        enable_full_relativity,
    )
    r_packet.nu = comoving_nu * inverse_doppler_factor
    r_packet.energy = comoving_energy * inverse_doppler_factor
    r_packet.mu = comoving_mu
    if enable_full_relativity:
        velocity = geometry.get_velocity(
            r_packet.r,
            r_packet.current_shell_id,
        )
        r_packet.mu = angle_aberration_CMF_to_LF_from_velocity(
            velocity,
            comoving_mu,
        )
