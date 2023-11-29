import logging
import os

import numpy as np
import pandas as pd
from astropy import units as u

from tardis.io.model.parse_density_configuration import (
    calculate_density_after_time,
    parse_config_v1_density,
    parse_csvy_density,
)
from tardis.io.model.readers.base import read_abundances_file, read_density_file
from tardis.io.model.readers.csvy import parse_csv_abundances
from tardis.io.model.readers.generic_readers import read_uniform_abundances
from tardis.model.geometry.radial1d import HomologousRadial1DGeometry
from tardis.model.matter.composition import Composition
from tardis.model.matter.decay import IsotopicMassFraction
from tardis.util.base import quantity_linspace

logger = logging.getLogger(__name__)


def parse_structure_config(config, time_explosion, enable_homology=True):
    """
    Parse the structure configuration data.

    Parameters
    ----------
    config : object
        The configuration data.
    time_explosion : float
        The time of the explosion.
    enable_homology : bool, optional
        Whether to enable homology (default is True).

    Returns
    -------
    electron_densities : object
        The parsed electron densities.
    temperature : object
        The parsed temperature.
    geometry : object
        The parsed geometry.
    density : object
        The parsed density.

    Raises
    ------
    NotImplementedError
        If the structure configuration type is not supported.

    Notes
    -----
    This function parses the structure configuration data and returns the parsed electron
    densities, temperature, geometry, and density. The structure configuration can be of
    type 'specific' or 'file'. If it is of type 'specific', the velocity and density are
    parsed from the configuration. If it is of type 'file', the velocity and density are
    read from a file. The parsed data is used to create a homologous radial 1D geometry object.
    """
    electron_densities = None
    temperature = None
    structure_config = config.model.structure
    if structure_config.type == "specific":
        velocity = quantity_linspace(
            structure_config.velocity.start,
            structure_config.velocity.stop,
            structure_config.velocity.num + 1,
        ).cgs
        density = parse_config_v1_density(config)

    elif structure_config.type == "file":
        if os.path.isabs(structure_config.filename):
            structure_config_fname = structure_config.filename
        else:
            structure_config_fname = os.path.join(
                config.config_dirname, structure_config.filename
            )

        (
            time_0,
            velocity,
            density_0,
            electron_densities,
            temperature,
        ) = read_density_file(structure_config_fname, structure_config.filetype)
        density_0 = density_0.insert(0, 0)

        density = calculate_density_after_time(
            density_0, time_0, time_explosion
        )

    else:
        raise NotImplementedError

    # Note: This is the number of shells *without* taking in mind the
    #       v boundaries.
    if len(density) == len(velocity):
        logger.warning(
            "Number of density points larger than number of shells. Assuming inner point irrelevant"
        )
        density = density[1:]
    geometry = HomologousRadial1DGeometry(
        velocity[:-1],  # r_inner
        velocity[1:],  # r_outer
        v_inner_boundary=structure_config.get("v_inner_boundary", None),
        v_outer_boundary=structure_config.get("v_outer_boundary", None),
        time_explosion=time_explosion,
    )
    return electron_densities, temperature, geometry, density


def parse_csvy_geometry(
    config, csvy_model_config, csvy_model_data, time_explosion
):
    """
    Parse the geometry data from a CSVY model.

    Parameters
    ----------
    config : object
        The configuration data.
    csvy_model_config : object
        The configuration data of the CSVY model.
    csvy_model_data : object
        The data of the CSVY model.
    time_explosion : float
        The time of the explosion.

    Returns
    -------
    geometry : object
        The parsed geometry.

    Raises
    ------
    None.

    Notes
    -----
    This function parses the geometry data from a CSVY model. It extracts the velocity
    information from the CSVY model configuration or data. The parsed velocity data is
    used to create a homologous radial 1D geometry object, which is returned.
    """
    if hasattr(config, "model"):
        if hasattr(config.model, "v_inner_boundary"):
            v_boundary_inner = config.model.v_inner_boundary
        else:
            v_boundary_inner = None

        if hasattr(config.model, "v_outer_boundary"):
            v_boundary_outer = config.model.v_outer_boundary
        else:
            v_boundary_outer = None
    else:
        v_boundary_inner = None
        v_boundary_outer = None

    if hasattr(csvy_model_config, "velocity"):
        velocity = quantity_linspace(
            csvy_model_config.velocity.start,
            csvy_model_config.velocity.stop,
            csvy_model_config.velocity.num + 1,
        ).cgs
    else:
        velocity_field_index = [
            field["name"] for field in csvy_model_config.datatype.fields
        ].index("velocity")
        velocity_unit = u.Unit(
            csvy_model_config.datatype.fields[velocity_field_index]["unit"]
        )
        velocity = csvy_model_data["velocity"].values * velocity_unit
        velocity = velocity.to("cm/s")

    geometry = HomologousRadial1DGeometry(
        velocity[:-1],  # r_inner
        velocity[1:],  # r_outer
        v_inner_boundary=v_boundary_inner,
        v_outer_boundary=v_boundary_outer,
        time_explosion=time_explosion,
    )
    return geometry


def parse_abundance_config(config, geometry, time_explosion):
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

    isotope_abundance = IsotopicMassFraction(
        isotope_abundance, time_0=config.model.abundances.model_isotope_time_0
    ).decay(time_explosion)

    nuclide_mass_fraction = convert_to_nuclide_mass_fraction(
        isotope_abundance, abundance
    )
    return nuclide_mass_fraction


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


def parse_composition_csvy(
    atom_data, csvy_model_config, csvy_model_data, time_explosion, geometry
):
    """
    Parse the composition data from a CSVY model.

    Parameters
    ----------
    atom_data : object
        The atom data used for parsing.
    csvy_model_config : object
        The configuration data of the CSVY model.
    csvy_model_data : object
        The data of the CSVY model.
    time_explosion : float
        The time of the explosion.
    geometry : object
        The geometry of the model.

    Returns
    -------
    density : object
        The parsed density data.
    abundance : object
        The parsed abundance data.
    isotope_abundance : object
        The parsed isotope abundance data.
    elemental_mass : object
        The elemental mass data.

    Raises
    ------
    None.

    Notes
    -----
    This function parses the composition data from a CSVY model. It calls the 'parse_density_csvy'
    function to parse the density data, and the 'parse_abundance_csvy' function to parse the abundance
    and isotope abundance data. The parsed data is returned as density, abundance, isotope_abundance,
    and elemental_mass.
    """
    density = parse_density_csvy(
        csvy_model_config, csvy_model_data, time_explosion
    )

    nuclide_mass_fraction = parse_abundance_csvy(
        csvy_model_config, csvy_model_data, geometry, time_explosion
    )
    return Composition(
        density, nuclide_mass_fraction, atom_data.atom_data.mass.copy()
    )


def parse_abundance_csvy(
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

    isotope_mass_fraction = IsotopicMassFraction(
        isotope_mass_fraction, time_0=csvy_model_config.model_isotope_time_0
    ).decay(time_explosion)
    return convert_to_nuclide_mass_fraction(
        isotope_mass_fraction, mass_fraction
    )


def parse_density_csvy(csvy_model_config, csvy_model_data, time_explosion):
    """
    Parse the density data from a CSVY model.

    Parameters
    ----------
    csvy_model_config : object
        The configuration data of the CSVY model.
    csvy_model_data : object
        The data of the CSVY model.
    time_explosion : float
        The time of the explosion.

    Returns
    -------
    density : object
        The parsed density data.

    Raises
    ------
    None.

    Notes
    -----
    This function parses the density data from a CSVY model. If the CSVY model configuration
    contains a 'density' attribute, it uses the 'parse_csvy_density' function to parse the
    density data. Otherwise, it calculates the density data using the 'calculate_density_after_time'
    function. The parsed density data is returned.
    """
    if hasattr(csvy_model_config, "density"):
        density = parse_csvy_density(csvy_model_config, time_explosion)
    else:
        time_0 = csvy_model_config.model_density_time_0
        density_field_index = [
            field["name"] for field in csvy_model_config.datatype.fields
        ].index("density")
        density_unit = u.Unit(
            csvy_model_config.datatype.fields[density_field_index]["unit"]
        )
        density_0 = csvy_model_data["density"].values * density_unit

        # Removing the 0-th shell from the density array in the CSVY case
        if density_0[0] >= 0 * u.g / u.cm**3:
            logger.warning(
                "TARDIS ignores the 0-th shell for the CSVY format. Make density of 0-th shell negative to remove this warning."
            )
        density_0 = density_0.to("g/cm^3")[1:]
        # density_0 = density_0.insert(0, 0)
        density = calculate_density_after_time(
            density_0, time_0, time_explosion
        )

    return density
