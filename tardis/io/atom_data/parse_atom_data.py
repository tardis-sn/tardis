import logging
from pathlib import Path

from tardis.io.atom_data.base import AtomData

logger = logging.getLogger(__name__)


def parse_atom_data(config, atom_data=None):
    """
    Parse atom data for the simulation.

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
                "please see https://github.com/tardis-sn/tardis-regression-data/tree/master/atom_data"
                " for the most recent version.",
            )
            raise

    return atom_data
