from tardis.io.model.readers.cmfgen_deprecated import (
    read_cmfgen_composition,
    read_cmfgen_density,
)
from tardis.io.model.readers.generic_readers import (
    ConfigurationError,
    read_csv_composition,
    read_simple_ascii_mass_fractions,
    read_simple_ascii_density,
)


import numpy as np
import pandas as pd

from tardis.io.model.readers.artis import read_artis_density


def read_mass_fractions_file(
    mass_fractions_filename,
    mass_fractions_filetype,
    inner_boundary_index=None,
    outer_boundary_index=None,
):
    """
    read different density file formats

    Parameters
    ----------
    mass_fractions_filename : str
        filename or path of the density file
    mass_fractions_filetype : str
        type of the density file
    inner_boundary_index : int
        index of the inner shell, default None
    outer_boundary_index : int
        index of the outer shell, default None
    """

    file_parsers = {
        "simple_ascii": read_simple_ascii_mass_fractions,
        "artis": read_simple_ascii_mass_fractions,
        "cmfgen_model": read_cmfgen_composition,
        "custom_composition": read_csv_composition,
    }

    isotope_mass_fractions = pd.DataFrame()
    if mass_fractions_filetype in ["cmfgen_model", "custom_composition"]:
        index, mass_fractions, isotope_mass_fractions = file_parsers[
            mass_fractions_filetype
        ](mass_fractions_filename)
    else:
        index, mass_fractions = file_parsers[mass_fractions_filetype](
            mass_fractions_filename
        )

    if outer_boundary_index is not None:
        outer_boundary_index_m1 = outer_boundary_index - 1
    else:
        outer_boundary_index_m1 = None
    index = index[inner_boundary_index:outer_boundary_index]
    mass_fractions = mass_fractions.loc[
        :, slice(inner_boundary_index, outer_boundary_index_m1)
    ]
    mass_fractions.columns = np.arange(len(mass_fractions.columns))
    return index, mass_fractions, isotope_mass_fractions


def read_density_file(filename, filetype):
    """
    read different density file formats

    Parameters
    ----------
    filename : str
        filename or path of the density file
    filetype : str
        type of the density file

    Returns
    -------
    time_of_model : astropy.units.Quantity
        time at which the model is valid
    velocity : np.ndarray
        the array containing the velocities
    unscaled_mean_densities : np.ndarray
        the array containing the densities
    electron_densities : np.ndarray
        The array containing electron densities
    temperature : np.ndarray
        The array containing temperatures
    """
    file_parsers = {
        "artis": read_artis_density,
        "simple_ascii": read_simple_ascii_density,
        "cmfgen_model": read_cmfgen_density,
    }

    electron_densities = None
    temperature = None
    if filetype == "cmfgen_model":
        (
            time_of_model,
            velocity,
            unscaled_mean_densities,
            electron_densities,
            temperature,
        ) = read_cmfgen_density(filename)
    else:
        (time_of_model, velocity, unscaled_mean_densities) = file_parsers[
            filetype
        ](filename)

    v_inner = velocity[:-1]
    v_outer = velocity[1:]

    invalid_volume_mask = (v_outer - v_inner) <= 0
    if invalid_volume_mask.sum() > 0:
        message = "\n".join(
            [
                f"cell {i:d}: v_inner {v_inner_i:s}, v_outer " f"{v_outer_i:s}"
                for i, v_inner_i, v_outer_i in zip(
                    np.arange(len(v_outer))[invalid_volume_mask],
                    v_inner[invalid_volume_mask],
                    v_outer[invalid_volume_mask],
                )
            ]
        )
        raise ConfigurationError(
            "Invalid volume of following cell(s):\n" f"{message:s}"
        )

    return (
        time_of_model,
        velocity,
        unscaled_mean_densities,
        electron_densities,
        temperature,
    )
