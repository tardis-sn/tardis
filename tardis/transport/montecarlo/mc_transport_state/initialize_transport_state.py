from astropy import units as u

from tardis.model.geometry.radial1d import HomologousRadial1DGeometry
from tardis.opacities.macro_atom.macroatom_state import MacroAtomState
from tardis.opacities.opacity_state import OpacityState
from tardis.plasma.base import BasePlasma
from tardis.transport.montecarlo.estimators.radfield_mc_estimators import (
    initialize_estimator_statistics,
)
from tardis.transport.montecarlo.mc_transport_state.base import (
    MonteCarloTransportState,
)
from tardis.transport.montecarlo.packet_source.base import BasePacketSource


def initialize_transport_state(
    geometry_state: HomologousRadial1DGeometry,
    time_explosion: u.Quantity,
    opacity_state: OpacityState,
    macro_atom_state: MacroAtomState,
    plasma: BasePlasma,
    packet_source: BasePacketSource,
    no_of_packets: int,
    line_interaction_type: str,
    iteration: int = 0,
) -> MonteCarloTransportState:
    """
    Initialize a MonteCarloTransportState from simulation components.

    Parameters
    ----------
    geometry_state : HomologousRadial1DGeometry
        Geometry state containing shell boundaries and velocity information.
    time_explosion : u.Quantity
        Time since explosion.
    opacity_state : OpacityState
        Opacity state containing line and continuum opacities.
    macro_atom_state : MacroAtomState
        Macro atom state with transition probabilities.
    plasma : BasePlasma
        Plasma state with physical properties.
    packet_source : BasePacketSource
        Packet source for creating initial packets.
    no_of_packets : int
        Number of packets to create.
    line_interaction_type : str
        Type of line interaction (e.g., 'scatter', 'macroatom', 'downbranch').
    iteration : int, optional
        Iteration number for seed offset. Default is 0.

    Returns
    -------
    MonteCarloTransportState
        Initialized transport state ready for Monte Carlo simulation.
    """
    if not plasma.continuum_interaction_species.empty:
        if plasma.gamma is not None:
            gamma_shape = plasma.gamma.shape
        else:
            gamma_shape = plasma.phi_lucy.shape

    else:
        gamma_shape = (0, 0)

    packet_collection = packet_source.create_packets(
        no_of_packets, seed_offset=iteration
    )

    geometry_state_numba = geometry_state.to_numba()
    opacity_state_numba = opacity_state.to_numba(
        macro_atom_state,
        line_interaction_type,
    )
    opacity_state_numba = opacity_state_numba[
        geometry_state.v_inner_boundary_index : geometry_state.v_outer_boundary_index
    ]

    estimators = initialize_estimator_statistics(
        opacity_state_numba.tau_sobolev.shape, gamma_shape
    )

    return MonteCarloTransportState(
        packet_collection,
        estimators,
        geometry_state=geometry_state_numba,
        opacity_state=opacity_state_numba,
        time_explosion=time_explosion,
    )
