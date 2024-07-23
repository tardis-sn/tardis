import logging
import os

import astropy.units as u
import numpy as np
import pandas as pd

from tardis.io.model.readers.base import read_mass_fractions_file
from tardis.io.model.readers.csvy import parse_csv_mass_fractions
from tardis.io.model.readers.generic_readers import read_uniform_mass_fractions
from tardis.model.matter.composition import Composition
from tardis.model.matter.decay import IsotopicMassFraction

logger = logging.getLogger(__name__)


def parse_mass_fractions_from_config(config, geometry, time_explosion):
    """
    Parse the mass fraction configuration data.

    Parameters
    ----------
    config : object
        The configuration data.
    geometry : object
        The geometry of the model.
    time_explosion : float
        The time of the explosion.

    Returns
    -------
    nuclide_mass_fractions : object
        The parsed nuclide mass fraction.

    raw_isotope_mass_fractions : object
        The parsed raw isotope mass fraction. This is the isotope mass fraction data before decay.

    Raises
    ------
    None.

    Notes
    -----
    This function parses the abundance configuration data and returns the parsed nuclide
    mass fraction. The mass fraction configuration can be of type 'uniform' or 'file'. If it
    is of type 'uniform', the mass fraction and isotope mass fraction are read using the
    'read_uniform_mass fractions' function. If it is of type 'file', the mass fraction and
    isotope mass fraction are read from a file using the 'read_abundances_file' function.
    The parsed data is then processed to replace NaN values with 0.0, remove rows with
    zero sum, and normalize the data if necessary. The resulting nuclide mass fraction
    is returned.
    """
    mass_fractions_section = config.model.abundances
    isotope_mass_fractions = pd.DataFrame()

    if mass_fractions_section.type == "uniform":
        mass_fractions, isotope_mass_fractions = read_uniform_mass_fractions(
            mass_fractions_section, geometry.no_of_shells
        )

    elif mass_fractions_section.type == "file":
        if os.path.isabs(mass_fractions_section.filename):
            mass_fractions_fname = mass_fractions_section.filename
        else:
            mass_fractions_fname = os.path.join(
                config.config_dirname, mass_fractions_section.filename
            )

        (
            index,
            mass_fractions,
            isotope_mass_fractions,
        ) = read_mass_fractions_file(
            mass_fractions_fname, mass_fractions_section.filetype
        )

    mass_fractions = mass_fractions.replace(np.nan, 0.0)
    mass_fractions = mass_fractions[mass_fractions.sum(axis=1) > 0]

    norm_factor = mass_fractions.sum(axis=0) + isotope_mass_fractions.sum(
        axis=0
    )

    if np.any(np.abs(norm_factor - 1) > 1e-12):
        logger.warning(
            "Mass fractions have not been normalized to 1. - normalizing"
        )
        mass_fractions /= norm_factor
        isotope_mass_fractions /= norm_factor
    # The next line is if the mass_fractions are given via dict
    # and not gone through the schema validator
    raw_isotope_mass_fractions = isotope_mass_fractions
    model_isotope_time_0 = config.model.abundances.get(
        "model_isotope_time_0", 0.0 * u.day
    )
    isotope_mass_fractions = IsotopicMassFraction(
        isotope_mass_fractions, time_0=model_isotope_time_0
    ).decay(time_explosion)

    nuclide_mass_fractions = convert_to_nuclide_mass_fractions(
        isotope_mass_fractions, mass_fractions
    )
    return nuclide_mass_fractions, raw_isotope_mass_fractions


def parse_mass_fractions_from_csvy(
    csvy_model_config, csvy_model_data, geometry, time_explosion
):
    """
    Parse the mass fraction data from a CSVY model.

    Parameters
    ----------
    csvy_model_config : object
        The configuration data of the CSVY model.
    csvy_model_data : object
        The data of the CSVY model.
    geometry : object
        The geometry of the model.

    Returns
    -------
    mass_fractions : pd.DataFrame
        The parsed mass fraction data.
    isotope_mass_fractions : pandas.DataFrame
        The parsed isotope mass fraction data.

    Raises
    ------
    None.

    Notes
    -----
    This function parses the mass fraction data from a CSVY model. If the CSVY model
    configuration contains an 'abundance' attribute, it uses the 'read_uniform_mass_fractions'
    function to parse the mass fraction and isotope mass fraction data. Otherwise, it uses the
    'parse_csv_mass_fractions' function to parse the data. The parsed data is then processed
    to replace NaN values with 0.0, remove rows with zero sum, and normalize the data
    if necessary. The resulting mass fraction and isotope mass fraction arrays are returned.
    """
    if hasattr(csvy_model_config, "abundance"):
        mass_fractions_section = csvy_model_config.abundance
        mass_fractions, isotope_mass_fractions = read_uniform_mass_fractions(
            mass_fractions_section, geometry.no_of_shells
        )
    else:
        _, mass_fractions, isotope_mass_fractions = parse_csv_mass_fractions(
            csvy_model_data
        )
        mass_fractions = mass_fractions.loc[:, 1:]
        mass_fractions.columns = np.arange(mass_fractions.shape[1])
        isotope_mass_fractions = isotope_mass_fractions.loc[:, 1:]
        isotope_mass_fractions.columns = np.arange(
            isotope_mass_fractions.shape[1]
        )

    mass_fractions = mass_fractions.replace(np.nan, 0.0)
    mass_fractions = mass_fractions[mass_fractions.sum(axis=1) > 0]
    isotope_mass_fractions = isotope_mass_fractions.replace(np.nan, 0.0)
    isotope_mass_fractions = isotope_mass_fractions[
        isotope_mass_fractions.sum(axis=1) > 0
    ]
    norm_factor = mass_fractions.sum(axis=0) + isotope_mass_fractions.sum(
        axis=0
    )

    if np.any(np.abs(norm_factor - 1) > 1e-12):
        logger.warning(
            "Mass fractions have not been normalized to 1. - normalizing"
        )
        mass_fractions /= norm_factor
        isotope_mass_fractions /= norm_factor

    raw_isotope_mass_fraction = isotope_mass_fractions
    isotope_mass_fractions = IsotopicMassFraction(
        isotope_mass_fractions, time_0=csvy_model_config.model_isotope_time_0
    ).decay(time_explosion)
    return (
        convert_to_nuclide_mass_fractions(
            isotope_mass_fractions, mass_fractions
        ),
        raw_isotope_mass_fraction,
    )


def convert_to_nuclide_mass_fractions(isotopic_mass_fractions, mass_fractions):
    """
    Convert the mass fraction and isotope mass fraction data to nuclide mass fraction.

    Parameters
    ----------
    isotope_mass_fraction : pandas.DataFrame
        The isotope mass fraction data.
    mass_fractions : pandas.DataFrame
        The mass fraction data.

    Returns
    -------
    nuclide_mass_fraction : pandas.DataFrame
        The converted nuclide mass fraction.

    Raises
    ------
    None.

    Notes
    -----
    This function converts the mass fraction and isotope mass fraction data to nuclide mass fraction.
    If the mass fraction data is not None, it is converted to nuclide mass fraction by mapping
    the mass fraction index to nuclide indices using the 'convert_element2nuclide_index' function.
    The resulting mass fraction data is then concatenated with the isotope mass fraction data to
    obtain the final nuclide mass fraction.
    """
    nuclide_mass_fractions = pd.DataFrame()
    if mass_fractions is not None:
        mass_fractions.index = Composition.convert_element2nuclide_index(
            mass_fractions.index
        )
        nuclide_mass_fractions = mass_fractions
    else:
        nuclide_mass_fractions = pd.DataFrame()

    if isotopic_mass_fractions is not None:
        nuclide_mass_fractions = pd.concat(
            [nuclide_mass_fractions, isotopic_mass_fractions]
        )
    return nuclide_mass_fractions
