import logging
from pathlib import Path

from tardis.io.atom_data.base import AtomData
from tardis.model import SimulationState
from tardis.model.parse_input import initialize_packet_source

logger = logging.getLogger(__name__)


def initialize_atom_data(config, atom_data=None):
    """
    Initialize atom data for the simulation.

    Parameters
    ----------
    config : object
        The configuration object containing information about the atom data.
    atom_data : object, optional
        Existing atom data to be used, if provided.

    Returns
    -------
    object
        The initialized atom data.

    Raises
    ------
    ValueError
        If no atom_data option is found in the configuration.
    """
    if atom_data is None:
        if "atom_data" in config:
            if Path(config.atom_data).is_absolute():
                atom_data_fname = Path(config.atom_data)
            else:
                atom_data_fname = Path(config.config_dirname) / config.atom_data

        else:
            raise ValueError("No atom_data option found in the configuration.")

        logger.info(f"\n\tReading Atomic Data from {atom_data_fname}")

        try:
            atom_data = AtomData.from_hdf(atom_data_fname)
        except TypeError as e:
            print(
                e,
                "Error might be from the use of an old-format of the atomic database, \n"
                "please see https://github.com/tardis-sn/tardis-refdata/tree/master/atom_data"
                " for the most recent version.",
            )
            raise

    return atom_data


def initialize_simulation_state(
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
