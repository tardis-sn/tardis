import logging
import os

import astropy.units as u
import numpy as np
import pandas as pd

from tardis.io.model.readers.base import read_abundances_file
from tardis.io.model.readers.csvy import parse_csv_abundances
from tardis.io.model.readers.generic_readers import read_uniform_abundances
from tardis.model.matter.composition import Composition
from tardis.model.matter.decay import IsotopicMassFraction

logger = logging.getLogger(__name__)


def parse_abundance_from_config(config, geometry, time_explosion):
    """
    Parse the abundance configuration data.

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
    nuclide_mass_fraction : object
        The parsed nuclide mass fraction.

    raw_isotope_abundance : object
        The parsed raw isotope abundance. This is the isotope abundance data before decay.

    Raises
    ------
    None.

    Notes
    -----
    This function parses the abundance configuration data and returns the parsed nuclide
    mass fraction. The abundance configuration can be of type 'uniform' or 'file'. If it
    is of type 'uniform', the abundance and isotope abundance are read using the
    'read_uniform_abundances' function. If it is of type 'file', the abundance and
    isotope abundance are read from a file using the 'read_abundances_file' function.
    The parsed data is then processed to replace NaN values with 0.0, remove rows with
    zero sum, and normalize the data if necessary. The resulting nuclide mass fraction
    is returned.
    """
    abundances_section = config.model.abundances
    isotope_abundance = pd.DataFrame()

    if abundances_section.type == "uniform":
        abundance, isotope_abundance = read_uniform_abundances(
            abundances_section, geometry.no_of_shells
        )

    elif abundances_section.type == "file":
        if os.path.isabs(abundances_section.filename):
            abundances_fname = abundances_section.filename
        else:
            abundances_fname = os.path.join(
                config.config_dirname, abundances_section.filename
            )

        index, abundance, isotope_abundance = read_abundances_file(
            abundances_fname, abundances_section.filetype
        )

    abundance = abundance.replace(np.nan, 0.0)
    abundance = abundance[abundance.sum(axis=1) > 0]

    norm_factor = abundance.sum(axis=0) + isotope_abundance.sum(axis=0)

    if np.any(np.abs(norm_factor - 1) > 1e-12):
        logger.warning(
            "Abundances have not been normalized to 1. - normalizing"
        )
        abundance /= norm_factor
        isotope_abundance /= norm_factor
    # The next line is if the abundances are given via dict
    # and not gone through the schema validator
    raw_isotope_abundance = isotope_abundance
    model_isotope_time_0 = config.model.abundances.get(
        "model_isotope_time_0", 0.0 * u.day
    )
    isotope_abundance = IsotopicMassFraction(
        isotope_abundance, time_0=model_isotope_time_0
    ).decay(time_explosion)

    nuclide_mass_fraction = convert_to_nuclide_mass_fraction(
        isotope_abundance, abundance
    )
    return nuclide_mass_fraction, raw_isotope_abundance


def parse_abundance_from_csvy(
    csvy_model_config, csvy_model_data, geometry, time_explosion
):
    """
    Parse the abundance data from a CSVY model.

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
    abundance : pd.DataFrame
        The parsed abundance data.
    isotope_abundance : pandas.DataFrame
        The parsed isotope abundance data.

    Raises
    ------
    None.

    Notes
    -----
    This function parses the abundance data from a CSVY model. If the CSVY model
    configuration contains an 'abundance' attribute, it uses the 'read_uniform_abundances'
    function to parse the abundance and isotope abundance data. Otherwise, it uses the
    'parse_csv_abundances' function to parse the data. The parsed data is then processed
    to replace NaN values with 0.0, remove rows with zero sum, and normalize the data
    if necessary. The resulting abundance and isotope abundance arrays are returned.
    """
    if hasattr(csvy_model_config, "abundance"):
        abundances_section = csvy_model_config.abundance
        mass_fraction, isotope_mass_fraction = read_uniform_abundances(
            abundances_section, geometry.no_of_shells
        )
    else:
        _, mass_fraction, isotope_mass_fraction = parse_csv_abundances(
            csvy_model_data
        )
        mass_fraction = mass_fraction.loc[:, 1:]
        mass_fraction.columns = np.arange(mass_fraction.shape[1])
        isotope_mass_fraction = isotope_mass_fraction.loc[:, 1:]
        isotope_mass_fraction.columns = np.arange(
            isotope_mass_fraction.shape[1]
        )

    mass_fraction = mass_fraction.replace(np.nan, 0.0)
    mass_fraction = mass_fraction[mass_fraction.sum(axis=1) > 0]
    isotope_mass_fraction = isotope_mass_fraction.replace(np.nan, 0.0)
    isotope_mass_fraction = isotope_mass_fraction[
        isotope_mass_fraction.sum(axis=1) > 0
    ]
    norm_factor = mass_fraction.sum(axis=0) + isotope_mass_fraction.sum(axis=0)

    if np.any(np.abs(norm_factor - 1) > 1e-12):
        logger.warning(
            "Abundances have not been normalized to 1. - normalizing"
        )
        mass_fraction /= norm_factor
        isotope_mass_fraction /= norm_factor

    raw_isotope_mass_fraction = isotope_mass_fraction
    isotope_mass_fraction = IsotopicMassFraction(
        isotope_mass_fraction, time_0=csvy_model_config.model_isotope_time_0
    ).decay(time_explosion)
    return (
        convert_to_nuclide_mass_fraction(isotope_mass_fraction, mass_fraction),
        raw_isotope_mass_fraction,
    )


def convert_to_nuclide_mass_fraction(isotopic_mass_fraction, mass_fraction):
    """
    Convert the abundance and isotope abundance data to nuclide mass fraction.

    Parameters
    ----------
    isotope_abundance : pandas.DataFrame
        The isotope abundance data.
    abundance : pandas.DataFrame
        The abundance data.

    Returns
    -------
    nuclide_mass_fraction : pandas.DataFrame
        The converted nuclide mass fraction.

    Raises
    ------
    None.

    Notes
    -----
    This function converts the abundance and isotope abundance data to nuclide mass fraction.
    If the abundance data is not None, it is converted to nuclide mass fraction by mapping
    the abundance index to nuclide indices using the 'convert_element2nuclide_index' function.
    The resulting abundance data is then concatenated with the isotope abundance data to
    obtain the final nuclide mass fraction.
    """
    nuclide_mass_fraction = pd.DataFrame()
    if mass_fraction is not None:
        mass_fraction.index = Composition.convert_element2nuclide_index(
            mass_fraction.index
        )
        nuclide_mass_fraction = mass_fraction
    else:
        nuclide_mass_fraction = pd.DataFrame()

    if isotopic_mass_fraction is not None:
        nuclide_mass_fraction = pd.concat(
            [nuclide_mass_fraction, isotopic_mass_fraction]
        )
    return nuclide_mass_fraction
