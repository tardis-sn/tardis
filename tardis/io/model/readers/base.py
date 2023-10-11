from tardis.io.model.parse_density_configuration import (
    calculate_density_after_time,
    parse_config_v1_density,
)
import os
from tardis.io.model.readers.cmfgen import (
    read_cmfgen_composition,
    read_cmfgen_density,
)
from tardis.io.model.readers.generic_readers import (
    ConfigurationError,
    read_csv_composition,
    read_simple_ascii_abundances,
    read_simple_ascii_density,
)


import numpy as np
import pandas as pd

from tardis.io.model.readers.artis import read_artis_density
from tardis.model.base import logger
from tardis.model.geometry.radial1d import HomologousRadial1DGeometry
from tardis.util.base import quantity_linspace


def read_abundances_file(
    abundance_filename,
    abundance_filetype,
    inner_boundary_index=None,
    outer_boundary_index=None,
):
    """
    read different density file formats

    Parameters
    ----------
    abundance_filename : str
        filename or path of the density file
    abundance_filetype : str
        type of the density file
    inner_boundary_index : int
        index of the inner shell, default None
    outer_boundary_index : int
        index of the outer shell, default None
    """

    file_parsers = {
        "simple_ascii": read_simple_ascii_abundances,
        "artis": read_simple_ascii_abundances,
        "cmfgen_model": read_cmfgen_composition,
        "custom_composition": read_csv_composition,
    }

    isotope_abundance = pd.DataFrame()
    if abundance_filetype in ["cmfgen_model", "custom_composition"]:
        index, abundances, isotope_abundance = file_parsers[abundance_filetype](
            abundance_filename
        )
    else:
        index, abundances = file_parsers[abundance_filetype](abundance_filename)

    if outer_boundary_index is not None:
        outer_boundary_index_m1 = outer_boundary_index - 1
    else:
        outer_boundary_index_m1 = None
    index = index[inner_boundary_index:outer_boundary_index]
    abundances = abundances.loc[
        :, slice(inner_boundary_index, outer_boundary_index_m1)
    ]
    abundances.columns = np.arange(len(abundances.columns))
    return index, abundances, isotope_abundance


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


def parse_structure_config(config, time_explosion, enable_homology=True):
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
