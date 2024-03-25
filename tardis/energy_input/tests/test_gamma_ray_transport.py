from pathlib import Path

import astropy.units as u
import numpy as np
import numpy.testing as npt
import pytest
import radioactivedecay as rd
from radioactivedecay import converters

from tardis.energy_input.gamma_ray_channel import (
    calculate_total_decays,
    create_inventories_dict,
    create_isotope_dicts,
)
from tardis.io.configuration import config_reader
from tardis.model import SimulationState


@pytest.fixture(scope="module")
def gamma_ray_config(example_configuration_dir: Path):
    """
    Parameters
    ----------
    example_configuration_dir: Path to the configuration directory.

    Returns
    -------
    Tardis configuration
    """
    yml_path = (
        example_configuration_dir
        / "tardis_configv1_density_exponential_nebular_multi_isotope.yml"
    )

    return config_reader.Configuration.from_yaml(yml_path)


@pytest.fixture(scope="module")
def gamma_ray_simulation_state(gamma_ray_config, atomic_dataset):
    """
    Parameters
    ----------
    gamma_ray_config:

    Returns
    -------
    Tardis model
    """
    gamma_ray_config.model.structure.velocity.start = 1.0 * u.km / u.s
    gamma_ray_config.model.structure.density.rho_0 = 5.0e2 * u.g / u.cm**3
    gamma_ray_config.supernova.time_explosion = 150 * u.d

    return SimulationState.from_config(
        gamma_ray_config, atom_data=atomic_dataset
    )


@pytest.fixture(scope="module")
def gamma_ray_model_state(gamma_ray_simulation_state):
    """
    Parameters
    ----------
    gamma_ray_simulation_state: Tardis simulation state

    Returns
    -------
    Tardis model state
    """

    raw_isotope_abundance = (
        gamma_ray_simulation_state.composition.raw_isotope_abundance
    )
    composition = gamma_ray_simulation_state.composition
    cell_masses = composition.calculate_cell_masses(
        gamma_ray_simulation_state.geometry.volume
    )

    return raw_isotope_abundance, cell_masses
