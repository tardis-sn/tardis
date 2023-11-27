import os
import pytest
import numpy as np
import pandas as pd
import numpy.testing as npt
from pathlib import Path
import radioactivedecay as rd
from radioactivedecay import converters

import tardis
from tardis.model import SimulationState
from tardis.io.configuration import config_reader
from tardis.io.util import yaml_load_file, YAMLLoader
from tardis.energy_input.energy_source import (
    get_nuclear_lines_database,
    get_all_isotopes,
    setup_input_energy,
)
from tardis.energy_input.gamma_ray_transport import (
    calculate_shell_masses,
    create_isotope_dicts,
    create_inventories_dict,
    calculate_total_decays,
    calculate_average_energies,
    decay_chain_energies,
    calculate_energy_per_mass,
)
import astropy.units as u
import astropy.constants as c


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
def simulation_setup(gamma_ray_config):
    """
    Parameters
    ----------
    gamma_ray_config:

    Returns
    -------
    Tardis model
    """
    gamma_ray_config.model.structure.velocity.start = 1.0 * u.km / u.s
    gamma_ray_config.model.structure.density.rho_0 = 5.0e2 * u.g / (u.cm**3)
    gamma_ray_config.supernova.time_explosion = 1.0 * (u.d)
    model = SimulationState.from_config(gamma_ray_config)
    return model


def test_calculate_shell_masses(simulation_setup):
    """Function to test calculation of shell masses.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    """
    model = simulation_setup
    volume = model.volume.to("cm^3")
    density = model.density.to("g/cm^3")
    desired = (volume * density).to(u.g).value

    shell_masses = calculate_shell_masses(model).value
    npt.assert_almost_equal(shell_masses, desired)


@pytest.mark.parametrize("nuclide_name", ["Ni-56", "Fe-52", "Cr-48"])
def test_activity(simulation_setup, nuclide_name):
    """
    Function to test the decay of an atom in radioactivedecay with an analytical solution.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    nuclide_name: Name of the nuclide.
    """
    nuclide = rd.Nuclide(nuclide_name)
    model = simulation_setup
    t_half = nuclide.half_life() * u.s
    decay_constant = np.log(2) / t_half
    time_delta = 1.0 * u.s
    shell_masses = calculate_shell_masses(model)
    raw_isotope_abundance = model.raw_isotope_abundance
    raw_isotope_abundance_mass = raw_isotope_abundance.apply(
        lambda x: x * shell_masses, axis=1
    )
    mass = raw_isotope_abundance_mass.loc[nuclide.Z, nuclide.A][0]
    iso_dict = create_isotope_dicts(raw_isotope_abundance, shell_masses)
    inv_dict = create_inventories_dict(iso_dict)

    total_decays = calculate_total_decays(inv_dict, time_delta)
    actual = total_decays[0][nuclide.Z, nuclide.A][nuclide_name]

    isotopic_mass = nuclide.atomic_mass * (u.g)
    number_of_moles = mass * (u.g) / isotopic_mass
    number_of_atoms = (number_of_moles * c.N_A).value
    N1 = number_of_atoms * np.exp(-decay_constant * time_delta)
    expected = number_of_atoms - N1.value

    npt.assert_allclose(actual, expected)


@pytest.mark.parametrize("nuclide_name", ["Ni-56", "Fe-52", "Cr-48"])
def test_activity_chain(simulation_setup, nuclide_name):
    """
    Function to test two atom decay chain in radioactivedecay with an analytical solution.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    nuclide_name: Name of the nuclide.
    """
    nuclide = rd.Nuclide(nuclide_name)
    model = simulation_setup
    t_half = nuclide.half_life()
    decay_constant = np.log(2) / t_half
    time_delta = 1.0 * (u.d).to(u.s)
    progeny = rd.Nuclide(nuclide.progeny()[0])
    decay_constant_progeny = np.log(2) / progeny.half_life()
    shell_masses = calculate_shell_masses(model)
    raw_isotope_abundance = model.raw_isotope_abundance
    raw_isotope_abundance_mass = raw_isotope_abundance.apply(
        lambda x: x * shell_masses, axis=1
    )
    mass = raw_isotope_abundance_mass.loc[nuclide.Z, nuclide.A][0]
    iso_dict = create_isotope_dicts(raw_isotope_abundance, shell_masses)
    inv_dict = create_inventories_dict(iso_dict)
    total_decays = calculate_total_decays(inv_dict, time_delta)
    actual_parent = total_decays[0][nuclide.Z, nuclide.A][nuclide_name]
    actual_progeny = total_decays[0][nuclide.Z, nuclide.A][nuclide.progeny()[0]]
    isotopic_mass = nuclide.atomic_mass
    number_of_moles = mass / isotopic_mass
    number_of_atoms = number_of_moles * converters.AVOGADRO
    decayed_parent = number_of_atoms * np.exp(-decay_constant * time_delta)
    expected_parent = number_of_atoms * (
        1 - np.exp(-decay_constant * time_delta)
    )

    npt.assert_almost_equal(expected_parent, actual_parent)
    # npt.assert_allclose(actual_progeny, expected_progeny, rtol=1e-3)


@pytest.mark.parametrize("nuclide_name", ["Ni-56", "Fe-52", "Cr-48"])
def test_isotope_dicts(simulation_setup, nuclide_name):
    """
    Function to test if the right names for the isotopes are present as dictionary keys.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    nuclide_name: Name of the nuclide.
    """
    model = simulation_setup
    nuclide = rd.Nuclide(nuclide_name)
    raw_isotope_abundance = model.raw_isotope_abundance
    shell_masses = calculate_shell_masses(model)
    iso_dict = create_isotope_dicts(raw_isotope_abundance, shell_masses)

    Z, A = nuclide.Z, nuclide.A

    for shell_number, isotope_dict in iso_dict.items():
        isotope_dict_key = isotope_dict[Z, A]
        assert nuclide_name.replace("-", "") in isotope_dict_key.keys()


@pytest.mark.parametrize("nuclide_name", ["Ni-56", "Fe-52", "Cr-48"])
def test_inventories_dict(simulation_setup, nuclide_name):
    """
    Function to test if the inventories dictionary is created correctly.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    nuclide_name: Name of the nuclide.
    """
    model = simulation_setup
    nuclide = rd.Nuclide(nuclide_name)
    raw_isotope_abundance = model.raw_isotope_abundance
    shell_masses = calculate_shell_masses(model)
    iso_dict = create_isotope_dicts(raw_isotope_abundance, shell_masses)
    inventories_dict = create_inventories_dict(iso_dict)

    Z, A = nuclide.Z, nuclide.A
    raw_isotope_abundance_mass = raw_isotope_abundance.apply(
        lambda x: x * shell_masses, axis=1
    )

    mass = raw_isotope_abundance_mass.loc[Z, A][0]
    inventory = rd.Inventory({nuclide.nuclide: mass}, "g")
    assert inventories_dict[0][Z, A] == inventory


def test_average_energies(simulation_setup, atomic_dataset):
    """
    Function to test if the energy from each isotope is there in the list.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    atomic_dataset: Tardis atomic and nuclear dataset.
    """

    model = simulation_setup
    raw_isotope_abundance = model.raw_isotope_abundance
    gamma_ray_lines = atomic_dataset.decay_radiation_data

    all_isotope_names = get_all_isotopes(raw_isotope_abundance)

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
def test_decay_energy_chain(simulation_setup, atomic_dataset, nuclide_name):
    """
    This function tests if the decay energy is calculated correctly for a decay chain.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    atomic_dataset: Tardis atomic and nuclear dataset.
    nuclide_name: Name of the nuclide.
    """
    model = simulation_setup
    nuclide = rd.Nuclide(nuclide_name)
    raw_isotope_abundance = model.raw_isotope_abundance
    shell_masses = calculate_shell_masses(model)
    iso_dict = create_isotope_dicts(raw_isotope_abundance, shell_masses)
    inventories_dict = create_inventories_dict(iso_dict)
    gamma_ray_lines = atomic_dataset.decay_radiation_data
    all_isotope_names = get_all_isotopes(raw_isotope_abundance)
    Z, A = nuclide.Z, nuclide.A

    total_decays = calculate_total_decays(inventories_dict, 1.0 * u.s)

    (
        average_energies,
        average_positron_energies,
        gamma_ray_line_dict,
    ) = calculate_average_energies(raw_isotope_abundance, gamma_ray_lines)

    decay_chain_energy = decay_chain_energies(
        average_energies,
        total_decays,
    )

    expected = (
        total_decays[0][Z, A][nuclide_name] * average_energies[nuclide_name]
    )
    actual = decay_chain_energy[0][Z, A][nuclide_name]

    npt.assert_almost_equal(expected, actual)


def test_energy_per_mass(simulation_setup, atomic_dataset):
    """
    This function tests if the energy per mass has the same dimensions as the raw_isotope_abundance.
    This means for each decay chain we are calculating the energy per mass, by summing the energy from each isotope.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    atomic_dataset: Tardis atomic and nuclear dataset.
    """
    model = simulation_setup
    raw_isotope_abundance = model.raw_isotope_abundance
    shell_masses = calculate_shell_masses(model)
    iso_dict = create_isotope_dicts(raw_isotope_abundance, shell_masses)
    inventories_dict = create_inventories_dict(iso_dict)
    total_decays = calculate_total_decays(inventories_dict, 1.0 * u.s)

    gamma_ray_lines = atomic_dataset.decay_radiation_data
    average_energies = calculate_average_energies(
        raw_isotope_abundance, gamma_ray_lines
    )
    decay_energy = decay_chain_energies(
        average_energies[0],
        total_decays,
    )
    energy_per_mass = calculate_energy_per_mass(
        decay_energy, raw_isotope_abundance, shell_masses
    )
    # If the shape is not same that means the code is not working with multiple isotopes
    assert (
        energy_per_mass[0].shape == (raw_isotope_abundance * shell_masses).shape
    )
