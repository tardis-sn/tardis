import pytest
import numpy as np
from pathlib import Path
import astropy.units as u
import numpy.testing as npt
import radioactivedecay as rd
import astropy.constants as const
from radioactivedecay import converters

from tardis.model import SimulationState
from tardis.io.configuration import config_reader
from tardis.energy_input.energy_source import (
    get_nuclear_lines_database,
)
from tardis.energy_input.gamma_ray_channel import (
    create_isotope_dicts,
    create_inventories_dict,
    calculate_total_decays,
    create_isotope_decay_df,
)


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
    gamma_ray_config: Tardis configuration
    atomic_dataset: Tardis atomic-nuclear dataset

    Returns
    -------
    Tardis simulation state
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
def test_isotope_dicts(gamma_ray_simulation_state, nuclide_name):
    """
    Function to test if the right names for the isotopes are present as dictionary keys.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    nuclide_name: Name of the nuclide.
    """
    nuclide = rd.Nuclide(nuclide_name)
    raw_isotope_abundance = (
        gamma_ray_simulation_state.composition.raw_isotope_abundance
    )
    composition = gamma_ray_simulation_state.composition
    cell_masses = composition.calculate_cell_masses(
        gamma_ray_simulation_state.geometry.volume
    )
    iso_dict = create_isotope_dicts(raw_isotope_abundance, cell_masses)

    for isotope_dict in iso_dict.values():
        assert nuclide_name.replace("-", "") in isotope_dict.keys()


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
    raw_isotope_abundance = (
        gamma_ray_simulation_state.composition.raw_isotope_abundance
    )
    composition = gamma_ray_simulation_state.composition
    cell_masses = composition.calculate_cell_masses(
        gamma_ray_simulation_state.geometry.volume
    )

    iso_dict = create_isotope_dicts(raw_isotope_abundance, cell_masses)
    inventories_dict = create_inventories_dict(iso_dict)

    Z, A = nuclide.Z, nuclide.A
    raw_isotope_abundance_mass = raw_isotope_abundance.apply(
        lambda x: x * cell_masses, axis=1
    )

    mass = raw_isotope_abundance_mass.loc[Z, A][0]
    isotope_inventory = rd.Inventory({nuclide.nuclide: mass}, "g")

    if nuclide_name in inventories_dict[0].contents:
        assert (
            inventories_dict[0].contents[nuclide_name]
            == isotope_inventory.contents[nuclide_name]
        )


@pytest.mark.parametrize("nuclide_name", ["Ni-56"])
def test_mass_energy_conservation(
    gamma_ray_simulation_state, atomic_dataset, nuclide_name
):
    """
    Function to test if the mass-energy conservation is satisfied.
    Parameters
    ----------
    simulation_setup: A simulation setup which returns a model.
    atomic_dataset: Tardis atomic-nuclear dataset
    nuclide_name: Name of the nuclide."""

    nuclide = rd.Nuclide(nuclide_name)
    raw_isotope_abundance = (
        gamma_ray_simulation_state.composition.raw_isotope_abundance
    )
    composition = gamma_ray_simulation_state.composition
    cell_masses = composition.calculate_cell_masses(
        gamma_ray_simulation_state.geometry.volume
    )
    gamma_ray_lines = atomic_dataset.decay_radiation_data
    iso_dict = create_isotope_dicts(raw_isotope_abundance, cell_masses)
    inventories_dict = create_inventories_dict(iso_dict)
    total_decays = calculate_total_decays(inventories_dict, 1 * u.d)
    isotope_decay_df = create_isotope_decay_df(total_decays, gamma_ray_lines)

    groupby_df = isotope_decay_df.groupby(level=["shell_number", "isotope"])

    ni_56 = (
        groupby_df.get_group((0, "Ni56"))["energy_per_channel"].sum()
        * (u.keV).to(u.MeV)
        * u.MeV
    )

    neutrino_energy = 0.41 * u.MeV

    total_energy_actual = ni_56 + neutrino_energy

    c2 = const.c.to("cm/s") ** 2

    # calculate mass of 56Ni
    ni56 = rd.Nuclide("Ni56")
    ni56_atomic_mass = ni56.atomic_mass * (u.u).to(u.g) * u.g

    # calculate mass of 56Co
    co56 = rd.Nuclide("co56")
    co56_atomic_mass = co56.atomic_mass * (u.u).to(u.g) * u.g

    Q = (ni56_atomic_mass - co56_atomic_mass) * c2 * u.erg.to(u.MeV)

    np.testing.assert_allclose(total_energy_actual.value, Q.value, rtol=0.01)
