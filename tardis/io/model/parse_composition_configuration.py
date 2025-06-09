import logging

from astropy import units as u

from tardis.io.configuration.config_reader import Configuration
from tardis.io.model.parse_density_configuration import (
    calculate_density_after_time,
    parse_density_from_csvy,
    parse_density_section_config,
)
from tardis.io.model.parse_geometry_configuration import (
    parse_structure_from_config,
)
from tardis.io.model.parse_mass_fraction_configuration import (
    parse_mass_fractions_from_config,
    parse_mass_fractions_from_csvy,
)
from tardis.model.matter.composition import Composition

logger = logging.getLogger(__name__)


def parse_density_from_config(config: Configuration) -> u.Quantity:
    """
    Parse the configuration file and produce a density at
    time_explosion.

    Parameters
    ----------
    config : tardis.io.config_reader.Configuration

    Returns
    -------
    density: u.Quantity

    """
    time_explosion = config.supernova.time_explosion.cgs
    (
        density_time,
        velocity,
        density,
        electron_densities,
        temperature,
    ) = parse_structure_from_config(config)

    if density is None:
        adjusted_velocity = velocity.insert(0, 0)
        v_middle = adjusted_velocity[1:] * 0.5 + adjusted_velocity[:-1] * 0.5
        d_conf = config.model.structure.density
        density, density_time = parse_density_section_config(
            d_conf, v_middle, time_explosion
        )

    density = calculate_density_after_time(
        density, density_time, time_explosion
    )
    # Note: This is the number of shells *without* taking in mind the
    #       v boundaries.
    if len(density) == len(velocity):
        logger.warning(
            "Number of density points larger than number of shells. Assuming inner point irrelevant"
        )
        density = density[1:]

    return density, electron_densities


def parse_composition_from_config(atom_data, config, time_explosion, geometry):
    """Parse the composition data from a CSVY model.

    Parameters
    ----------
    atom_data : object
        The atom data used for parsing.
    config : object
        The configuration data.
    time_explosion : float
        The time of the explosion.
    geometry : object
        The geometry of the model.

    Returns
    -------
    Composition
        The parsed composition
    array
        Electron densities.
    """
    density, electron_densities = parse_density_from_config(config)

    nuclide_mass_fractions = parse_mass_fractions_from_config(
        config, geometry, time_explosion
    )

    return (
        Composition(
            density,
            nuclide_mass_fractions,
            atom_data.atom_data.mass.copy(),
        ),
        electron_densities,
    )


def parse_composition_from_csvy(
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
    Composition : object
        The parsed composition

    Raises
    ------
    None.

    Notes
    -----
    This function parses the composition data from a CSVY model. It calls the 'parse_density_csvy'
    function to parse the density data, and the 'parse_mass_fraction_csvy' function to parse the mass fraction
    and isotope mass fraction data. The parsed data is returned as a Composition object.
    """
    density = parse_density_from_csvy(
        csvy_model_config, csvy_model_data, time_explosion
    )

    nuclide_mass_fractions = parse_mass_fractions_from_csvy(
        csvy_model_config, csvy_model_data, geometry, time_explosion
    )
    return Composition(
        density,
        nuclide_mass_fractions,
        atom_data.atom_data.mass.copy(),
    )
