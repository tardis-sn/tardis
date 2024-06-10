from tardis.io.model.parse_abundance_configuration import parse_abundance_csvy
from tardis.io.model.parse_density_configuration import parse_density_csvy
from tardis.model.matter.composition import Composition


def parse_csvy_composition(
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

    nuclide_mass_fraction, raw_isotope_mass_fraction = parse_abundance_csvy(
        csvy_model_config, csvy_model_data, geometry, time_explosion
    )
    return Composition(
        density,
        nuclide_mass_fraction,
        raw_isotope_mass_fraction,
        atom_data.atom_data.mass.copy(),
    )
