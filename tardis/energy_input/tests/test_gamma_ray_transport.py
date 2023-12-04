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


def test_calculate_cell_masses(gamma_ray_simulation_state):
    """Function to test calculation of shell masses.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    """
    volume = 2.70936170e39 * u.cm**3
    density = 5.24801665e-09 * u.g / u.cm**3
    desired = volume * density

    shell_masses = gamma_ray_simulation_state.composition.calculate_cell_masses(
        gamma_ray_simulation_state.geometry.volume
    )

    npt.assert_allclose(shell_masses[0], desired)


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


@pytest.mark.parametrize("nuclide_name", ["Ni-56", "Fe-52", "Cr-48"])
def test_isotope_dicts(gamma_ray_simulation_state, nuclide_name):
    """
    Function to test if the right names for the isotopes are present as dictionary keys.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    nuclide_name: Name of the nuclide.
    """
    nuclide = rd.Nuclide(nuclide_name)
    isotopic_mass_fractions = (
        gamma_ray_simulation_state.composition.isotopic_mass_fraction
    )
    composition = gamma_ray_simulation_state.composition
    cell_masses = composition.calculate_cell_masses(
        gamma_ray_simulation_state.geometry.volume
    )
    iso_dict = create_isotope_dicts(isotopic_mass_fractions, cell_masses)

    Z, A = nuclide.Z, nuclide.A

    for isotope_dict in iso_dict.values():
        isotope_dict_key = isotope_dict[Z, A]
        assert nuclide_name.replace("-", "") in isotope_dict_key.keys()


@pytest.mark.parametrize("nuclide_name", ["Ni-56", "Fe-52", "Cr-48"])
def test_inventories_dict(gamma_ray_simulation_state, nuclide_name):
    """
    Function to test if the inventories dictionary is created correctly.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    nuclide_name: Name of the nuclide.
    """

    nuclide = rd.Nuclide(nuclide_name)
    isotopic_mass_fractions = (
        gamma_ray_simulation_state.composition.isotopic_mass_fraction
    )
    composition = gamma_ray_simulation_state.composition
    cell_masses = composition.calculate_cell_masses(
        gamma_ray_simulation_state.geometry.volume
    )

    iso_dict = create_isotope_dicts(isotopic_mass_fractions, cell_masses)
    inventories_dict = create_inventories_dict(iso_dict)

    Z, A = nuclide.Z, nuclide.A
    raw_isotope_abundance_mass = isotopic_mass_fractions.apply(
        lambda x: x * cell_masses, axis=1
    )

    mass = raw_isotope_abundance_mass.loc[Z, A][0]
    inventory = rd.Inventory({nuclide.nuclide: mass}, "g")
    assert inventories_dict[0][Z, A] == inventory


def test_average_energies(gamma_ray_simulation_state, atomic_dataset):
    """
    Function to test if the energy from each isotope is there in the list.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    atomic_dataset: Tardis atomic and nuclear dataset.
    """

    isotopic_mass_fraction = (
        gamma_ray_simulation_state.composition.isotopic_mass_fraction
    )
    gamma_ray_lines = atomic_dataset.decay_radiation_data

    all_isotope_names = get_all_isotopes(isotopic_mass_fraction)

    average_energies_list = []

    for isotope_name in all_isotope_names:
        energy, intensity = setup_input_energy(
            gamma_ray_lines[
                gamma_ray_lines.index == isotope_name.replace("-", "")
            ],
            "g",
        )
        average_energies_list.append(np.sum(energy * intensity))  # keV

    assert len(average_energies_list) == len(all_isotope_names)


@pytest.mark.parametrize("nuclide_name", ["Ni-56", "Fe-52", "Cr-48"])
def test_decay_energy_chain(
    gamma_ray_simulation_state, atomic_dataset, nuclide_name
):
    """
    This function tests if the decay energy is calculated correctly for a decay chain.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    atomic_dataset: Tardis atomic and nuclear dataset.
    nuclide_name: Name of the nuclide.
    """

    nuclide = rd.Nuclide(nuclide_name)
    isotopic_mass_fractions = (
        gamma_ray_simulation_state.composition.isotopic_mass_fraction
    )

    composition = gamma_ray_simulation_state.composition
    cell_masses = composition.calculate_cell_masses(
        gamma_ray_simulation_state.geometry.volume
    )
    iso_dict = create_isotope_dicts(isotopic_mass_fractions, cell_masses)
    inventories_dict = create_inventories_dict(iso_dict)
    gamma_ray_lines = atomic_dataset.decay_radiation_data

    Z, A = nuclide.Z, nuclide.A

    total_decays = calculate_total_decays(inventories_dict, 1.0 * u.s)

    (
        average_energies,
        _,
        _,
    ) = calculate_average_energies(isotopic_mass_fractions, gamma_ray_lines)

    decay_chain_energy = decay_chain_energies(
        average_energies,
        total_decays,
    )

    expected = (
        total_decays[0][Z, A][nuclide_name] * average_energies[nuclide_name]
    )
    actual = decay_chain_energy[0][Z, A][nuclide_name]

    npt.assert_almost_equal(expected, actual)


def test_energy_per_mass(gamma_ray_simulation_state, atomic_dataset):
    """
    This function tests if the energy per mass has the same dimensions as the raw_isotope_abundance.
    This means for each decay chain we are calculating the energy per mass, by summing the energy from each isotope.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    atomic_dataset: Tardis atomic and nuclear dataset.
    """

    isotopic_mass_fractions = (
        gamma_ray_simulation_state.composition.isotopic_mass_fraction
    )
    composition = gamma_ray_simulation_state.composition
    cell_masses = composition.calculate_cell_masses(
        gamma_ray_simulation_state.geometry.volume
    )
    iso_dict = create_isotope_dicts(isotopic_mass_fractions, cell_masses)
    inventories_dict = create_inventories_dict(iso_dict)
    total_decays = calculate_total_decays(inventories_dict, 1.0 * u.s)

    gamma_ray_lines = atomic_dataset.decay_radiation_data
    average_energies = calculate_average_energies(
        isotopic_mass_fractions, gamma_ray_lines
    )
    decay_energy = decay_chain_energies(
        average_energies[0],
        total_decays,
    )
    energy_per_mass = calculate_energy_per_mass(
        decay_energy, isotopic_mass_fractions, cell_masses
    )
    # If the shape is not same that means the code is not working with multiple isotopes
    assert (
        energy_per_mass[0].shape
        == (isotopic_mass_fractions * cell_masses).shape
    )
