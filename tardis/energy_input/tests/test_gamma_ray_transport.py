from pathlib import Path

import astropy.units as u
import numpy as np
import numpy.testing as npt
import pytest
import radioactivedecay as rd
from radioactivedecay import converters

from tardis.energy_input.energy_source import (
    get_all_isotopes,
    setup_input_energy,
)
from tardis.energy_input.gamma_ray_transport import (
    calculate_average_energies,
    calculate_energy_per_mass,
    calculate_total_decays,
    create_inventories_dict,
    create_isotope_dicts,
    decay_chain_energies,
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


@pytest.mark.parametrize("nuclide_name", ["Ni-56", "Fe-52", "Cr-48"])
def test_activity(gamma_ray_simulation_state, nuclide_name):
    """
    Function to test the decay of an atom in radioactivedecay with an analytical solution.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    nuclide_name: Name of the nuclide.
    """
    # setup of decay test
    nuclide = rd.Nuclide(nuclide_name)
    t_half = nuclide.half_life() * u.s
    decay_constant = np.log(2) / t_half
    time_delta = 1.0 * u.s

    # calculating necessary values
    composition = gamma_ray_simulation_state.composition
    cell_masses = composition.calculate_cell_masses(
        gamma_ray_simulation_state.geometry.volume
    )
    isotopic_mass_fractions = (
        gamma_ray_simulation_state.composition.isotopic_mass_fraction
    )
    isotopic_masses = isotopic_mass_fractions * cell_masses
    test_mass = isotopic_masses.loc[(nuclide.Z, nuclide.A), 0] * u.g
    iso_dict = create_isotope_dicts(isotopic_mass_fractions, cell_masses)
    inv_dict = create_inventories_dict(iso_dict)

    total_decays = calculate_total_decays(inv_dict, time_delta)
    actual = total_decays[0][nuclide.Z, nuclide.A][nuclide_name]

    isotope_mass = nuclide.atomic_mass * u.u
    number_of_atoms = (test_mass / isotope_mass).to(u.dimensionless_unscaled)
    expected = number_of_atoms * (1 - np.exp(-decay_constant * time_delta))

    npt.assert_allclose(actual, expected)


@pytest.mark.parametrize("nuclide_name", ["Ni-56", "Fe-52", "Cr-48"])
def test_activity_chain(gamma_ray_simulation_state, nuclide_name):
    """
    Function to test two atom decay chain in radioactivedecay with an analytical solution.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    nuclide_name: Name of the nuclide.
    """
    nuclide = rd.Nuclide(nuclide_name)
    t_half = nuclide.half_life()
    decay_constant = np.log(2) / t_half
    time_delta = 1.0 * (u.d).to(u.s)

    composition = gamma_ray_simulation_state.composition
    cell_masses = composition.calculate_cell_masses(
        gamma_ray_simulation_state.geometry.volume
    )
    isotopic_mass_fractions = (
        gamma_ray_simulation_state.composition.isotopic_mass_fraction
    )
    isotopic_masses = isotopic_mass_fractions * cell_masses
    test_mass = isotopic_masses.loc[(nuclide.Z, nuclide.A), 0] * u.g
    iso_dict = create_isotope_dicts(isotopic_mass_fractions, cell_masses)
    inv_dict = create_inventories_dict(iso_dict)

    total_decays = calculate_total_decays(inv_dict, time_delta)
    actual_parent = total_decays[0][nuclide.Z, nuclide.A][nuclide_name]

    isotopic_mass = nuclide.atomic_mass * u.g
    number_of_moles = test_mass / isotopic_mass
    number_of_atoms = number_of_moles * converters.AVOGADRO
    expected_parent = number_of_atoms.to(1).value * (
        1 - np.exp(-decay_constant * time_delta)
    )

    npt.assert_almost_equal(expected_parent, actual_parent)
