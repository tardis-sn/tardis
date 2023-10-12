import logging
import os

from astropy import units as u
from tardis.model.matter.decay import IsotopeAbundances
import numpy as np
import pandas as pd
from tardis.io.model.parse_density_configuration import (
    calculate_density_after_time,
    parse_config_v1_density,
)
from tardis.io.model.readers.base import read_abundances_file, read_density_file
from tardis.io.model.readers.generic_readers import read_uniform_abundances
from tardis.model.geometry.radial1d import HomologousRadial1DGeometry
from tardis.util.base import quantity_linspace

logger = logging.getLogger(__name__)


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


def parse_csvy_geometry(
    config, csvy_model_config, csvy_model_data, time_explosion
):
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


def parse_abundance_section(config, atom_data, geometry):
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
            "Abundances have not been normalized to 1." " - normalizing"
        )
        abundance /= norm_factor
        isotope_abundance /= norm_factor

    isotope_abundance = IsotopeAbundances(isotope_abundance)

    elemental_mass = None
    if atom_data is not None:
        elemental_mass = atom_data.atom_data.mass

    return isotope_abundance, abundance, elemental_mass
