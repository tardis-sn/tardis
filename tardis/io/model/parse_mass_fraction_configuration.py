import logging
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd

from tardis.io.configuration.config_reader import (
    Configuration,
)
from tardis.io.model.csvy import parse_csv_mass_fractions
from tardis.io.model.readers.base import read_mass_fractions_file
from tardis.io.model.readers.generic_readers import read_uniform_mass_fractions
from tardis.model.geometry.radial1d import HomologousRadial1DGeometry
from tardis.model.matter.composition import Composition
from tardis.model.matter.decay import IsotopicMassFraction

logger = logging.getLogger(__name__)


def parse_mass_fractions_from_config(
    config: Configuration,
    geometry: HomologousRadial1DGeometry,
    time_explosion: u.Quantity,
) -> pd.DataFrame:
    """Parse the mass fraction configuration data.

    Parameters
    ----------
    config
        The configuration data.
    geometry
        The geometry of the model.
    time_explosion
        The time of the explosion.

    Returns
    -------
    nuclide_mass_fractions
        The parsed nuclide mass fraction.

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
        if Path(mass_fractions_section.filename).is_absolute():
            mass_fractions_fname = mass_fractions_section.filename
        else:
            mass_fractions_fname = (
                Path(config.config_dirname) / mass_fractions_section.filename
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

    model_isotope_time_0 = config.model.abundances.model_isotope_time_0

    if not np.isnan(model_isotope_time_0):
        assert model_isotope_time_0 < time_explosion
        isotope_mass_fractions = IsotopicMassFraction(
            isotope_mass_fractions, time_0=model_isotope_time_0
        ).calculate_decayed_mass_fractions(time_explosion)
    else:
        logger.warning(
            "model_isotope_time_0 is not set in the configuration. "
            "Isotopic mass fractions will not be decayed and is assumed to be correct for the time_explosion. THIS IS NOT RECOMMENDED!"
        )
        isotope_mass_fractions = IsotopicMassFraction(isotope_mass_fractions)
    nuclide_mass_fractions = convert_to_nuclide_mass_fractions(
        isotope_mass_fractions, mass_fractions
    )
    return nuclide_mass_fractions


def parse_mass_fractions_from_csvy(
    csvy_model_config: Configuration,
    csvy_model_data: pd.DataFrame | None,
    geometry: HomologousRadial1DGeometry,
    time_explosion: u.Quantity,
) -> pd.DataFrame:
    """Parse the mass fraction data from a CSVY model.

    Parameters
    ----------
    csvy_model_config
        The configuration data of the CSVY model.
    csvy_model_data
        The data of the CSVY model.
    geometry
        The geometry of the model.
    time_explosion
        The time of the explosion.

    Returns
    -------
    nuclide_mass_fractions
        The parsed nuclide mass fraction data.

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

    isotope_mass_fractions = IsotopicMassFraction(
        isotope_mass_fractions, time_0=csvy_model_config.model_isotope_time_0
    ).calculate_decayed_mass_fractions(time_explosion)
    return convert_to_nuclide_mass_fractions(
        isotope_mass_fractions, mass_fractions
    )


def convert_to_nuclide_mass_fractions(
    isotopic_mass_fractions: pd.DataFrame | None,
    mass_fractions: pd.DataFrame | None,
) -> pd.DataFrame:
    """Convert the mass fraction and isotope mass fraction data to nuclide mass fraction.

    Parameters
    ----------
    isotopic_mass_fractions
        The isotope mass fraction data.
    mass_fractions
        The mass fraction data.

    Returns
    -------
    nuclide_mass_fractions
        The converted nuclide mass fraction.

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
