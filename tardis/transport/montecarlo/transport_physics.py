"""Resolve compatible transport physics before compiled packet transport."""

from collections.abc import Callable
from dataclasses import dataclass

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.model.geometry.radial1d_homologous import (
    NumbaHomologousRadial1DGeometry,
)
from tardis.transport.montecarlo.modes.classic.packet_propagation import (
    packet_propagation as propagate_packet_homologous,
)
from tardis.transport.montecarlo.modes.iip.montecarlo_transport import (
    montecarlo_transport as transport_packets_continuum,
)
from tardis.transport.montecarlo.modes.iip.packet_propagation import (
    packet_propagation as propagate_packet_continuum,
)
from tardis.transport.montecarlo.modes.montecarlo_transport import (
    montecarlo_transport_with_vpackets as transport_packets_with_vpackets,
)
from tardis.transport.montecarlo.modes.nonhomologous.packet_propagation import (
    packet_propagation as propagate_packet_piecewise_linear,
)

type GeometryState = NumbaHomologousRadial1DGeometry | NumbaRadial1DGeometry
type PacketPropagationFunction = Callable[..., None]
type TransportFunction = Callable[..., tuple]


@dataclass(frozen=True, slots=True)
class ResolvedTransportPhysics:
    """Store a compatible compiled transport kernel and frame treatment.

    Packet propagation is the temporary kinematic boundary because the current
    nonhomologous implementation changes line initialization, resonance
    ordering, cursor updates, estimators, virtual packets, and interaction
    events together. Once those operations share a common contract, this field
    can become a smaller resonance callback without changing the resolver.
    ``transport_packets`` and ``propagate_packet`` are a matched pair and must
    not be recombined independently.

    Attributes
    ----------
    transport_packets : collections.abc.Callable
        Compiled outer driver compatible with the resolved physics.
    propagate_packet : collections.abc.Callable
        Compiled packet state machine compatible with the resolved physics.
    enable_full_relativity : bool
        Whether to use TARDIS's existing full-relativity branch.
    """

    transport_packets: TransportFunction
    propagate_packet: PacketPropagationFunction
    enable_full_relativity: bool


def resolve_transport_physics(
    geometry: GeometryState,
    enable_full_relativity: bool,
    continuum_process_enabled: bool,
) -> ResolvedTransportPhysics:
    """Resolve a supported transport kernel before entering Numba.

    Parameters
    ----------
    geometry : NumbaHomologousRadial1DGeometry or NumbaRadial1DGeometry
        Geometry state defining the velocity field.
    enable_full_relativity : bool
        Whether to use TARDIS's existing full-relativity branch.
    continuum_process_enabled : bool
        Whether bound-free and free-free interactions are enabled.

    Returns
    -------
    ResolvedTransportPhysics
        Compatible compiled packet kernel and frame-treatment choice.

    Raises
    ------
    NotImplementedError
        If the requested physics combination is not implemented.
    TypeError
        If the geometry is not supported by Monte Carlo transport.
    """
    if isinstance(geometry, NumbaHomologousRadial1DGeometry):
        if continuum_process_enabled:
            if not enable_full_relativity:
                raise NotImplementedError(
                    "Continuum transport currently requires full relativity."
                )
            propagate_packet = propagate_packet_continuum
            transport_packets = transport_packets_continuum
        else:
            propagate_packet = propagate_packet_homologous
            transport_packets = transport_packets_with_vpackets
    elif isinstance(geometry, NumbaRadial1DGeometry):
        if enable_full_relativity:
            raise NotImplementedError(
                "Full relativity is not implemented for piecewise-linear "
                "kinematics."
            )
        if continuum_process_enabled:
            raise NotImplementedError(
                "Continuum transport is not implemented for piecewise-linear "
                "kinematics."
            )
        propagate_packet = propagate_packet_piecewise_linear
        transport_packets = transport_packets_with_vpackets
    else:
        raise TypeError(
            f"Unsupported transport geometry: {type(geometry).__name__}"
        )

    return ResolvedTransportPhysics(
        transport_packets,
        propagate_packet,
        enable_full_relativity,
    )
