from tardis.io.configuration.config_reader import Configuration
from tardis.io.model.parse_packet_source_configuration import (
    initialize_packet_source,
)
from tardis.model import SimulationState
from tardis.transport.montecarlo.packet_source.base import BasePacketSource


def parse_simulation_state(
    config: Configuration,
    packet_source: BasePacketSource | None,
    enable_legacy_mode: bool,
    kwargs: dict,
    atom_data,
) -> SimulationState:
    """Initialize the simulation state.

    Parameters
    ----------
    config
        The configuration object for the simulation.
    packet_source
        The packet source for the simulation.
    enable_legacy_mode
        Flag indicating if legacy mode is enabled.
    kwargs
        Additional keyword arguments.
    atom_data
        The atom data for the simulation.

    Returns
    -------
    simulation_state
        The initialized simulation state.
    """
    if "model" in kwargs:
        simulation_state = kwargs["model"]
    else:
        if hasattr(config, "csvy_model"):
            simulation_state = SimulationState.from_csvy(
                config, legacy_mode_enabled=enable_legacy_mode
            )
        else:
            simulation_state = SimulationState.from_config(
                config,
                atom_data=atom_data,
                legacy_mode_enabled=enable_legacy_mode,
            )
        if packet_source is not None:
            simulation_state.packet_source = initialize_packet_source(
                packet_source, config, simulation_state.geometry
            )

    return simulation_state
