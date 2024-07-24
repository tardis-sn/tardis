from tardis.io.model.parse_packet_source_configuration import (
    initialize_packet_source,
)
from tardis.model import SimulationState


def parse_simulation_state(
    config, packet_source, legacy_mode_enabled, kwargs, atom_data
):
    """
    Initialize the simulation state.

    Parameters
    ----------
    config : object
        The configuration object for the simulation.
    packet_source : object
        The packet source for the simulation.
    legacy_mode_enabled : bool
        Flag indicating if legacy mode is enabled.
    kwargs : dict
        Additional keyword arguments.
    atom_data : object
        The atom data for the simulation.

    Returns
    -------
    object
        The initialized simulation state.
    """
    if "model" in kwargs:
        simulation_state = kwargs["model"]
    else:
        if hasattr(config, "csvy_model"):
            simulation_state = SimulationState.from_csvy(
                config,
                atom_data=atom_data,
                legacy_mode_enabled=legacy_mode_enabled,
            )
        else:
            simulation_state = SimulationState.from_config(
                config,
                atom_data=atom_data,
                legacy_mode_enabled=legacy_mode_enabled,
            )
        if packet_source is not None:
            simulation_state.packet_source = initialize_packet_source(
                config,
                simulation_state.geometry,
                packet_source,
                legacy_mode_enabled,
            )

    return simulation_state
