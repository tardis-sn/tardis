from pathlib import Path

import astropy.constants as const
import astropy.units as u
import numpy as np
import numpy.testing as npt
import pytest
import radioactivedecay as rd

from tardis.energy_input.gamma_ray_channel import (
    calculate_total_decays,
    create_inventories_dict,
    create_isotope_decay_df,
    create_isotope_dicts,
    time_evolve_cumulative_decay,
)
from tardis.energy_input.gamma_ray_transport import get_taus
from tardis.energy_input.main_gamma_ray_loop import get_effective_time_array
from tardis.energy_input.util import KEV2ERG
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


@pytest.fixture(scope="module")
def gamma_ray_test_composition(gamma_ray_simulation_state):
    """
    Parameters
    ----------
    gamma_ray_simulation_state: Tardis simulation state

    Returns
    -------
    raw_isotopic_mass_fraction: Raw isotopic mass fraction
    cell_masses: Mass of the cell
    """

    raw_isotopic_mass_fraction = (
        gamma_ray_simulation_state.composition.raw_isotope_abundance
    )
    composition = gamma_ray_simulation_state.composition
    cell_masses = composition.calculate_cell_masses(
        gamma_ray_simulation_state.geometry.volume
    )

    return raw_isotopic_mass_fraction, cell_masses


def test_calculate_cell_masses(gamma_ray_simulation_state):
    """Function to test calculation of shell masses.
    Parameters
    ----------
    gamma_ray_simulation_state: Tardis simulation state.
    """
    volume = 2.70936170e39 * u.cm**3
    density = 5.24801665e-09 * u.g / u.cm**3
    desired = volume * density

    shell_masses = gamma_ray_simulation_state.composition.calculate_cell_masses(
        gamma_ray_simulation_state.geometry.volume
    )

    npt.assert_allclose(shell_masses[0], desired)


@pytest.mark.parametrize("nuclide_name", ["Ni56", "Fe52", "Cr48"])
def test_isotope_dicts(gamma_ray_test_composition, nuclide_name):
    """
    Function to test if the right names for the isotopes are present as dictionary keys.
    Parameters
    ----------
    gamma_ray_test_composition: Function holding the composition.
    nuclide_name: Name of the nuclide.
    """
    raw_isotopic_mass_fraction, cell_masses = gamma_ray_test_composition
    isotope_dict = create_isotope_dicts(raw_isotopic_mass_fraction, cell_masses)

    for isotope_dict in isotope_dict.values():
        assert nuclide_name in isotope_dict.keys()


@pytest.mark.parametrize("nuclide_name", ["Ni-56", "Fe-52", "Cr-48"])
def test_inventories_dict(gamma_ray_test_composition, nuclide_name):
    """
    Function to test if the inventories dictionary is created correctly.
    Parameters
    ----------
    gamma_ray_test_composition: Function holding the composition.
    nuclide_name: Name of the nuclide.
    """

    nuclide = rd.Nuclide(nuclide_name)
    raw_isotopic_mass_fraction, cell_masses = gamma_ray_test_composition
    isotope_dict = create_isotope_dicts(raw_isotopic_mass_fraction, cell_masses)
    inventories_dict = create_inventories_dict(isotope_dict)

    Z, A = nuclide.Z, nuclide.A
    raw_isotope_mass = raw_isotopic_mass_fraction.apply(
        lambda x: x * cell_masses, axis=1
    )

    mass = raw_isotope_mass.loc[Z, A][0]
    isotope_inventory = rd.Inventory({nuclide.nuclide: mass}, "g")

    if nuclide_name in inventories_dict[0].contents:
        assert (
            inventories_dict[0].contents[nuclide_name]
            == isotope_inventory.contents[nuclide_name]
        )


@pytest.mark.parametrize("nuclide_name", ["Ni-56"])
def test_mass_energy_conservation(
    gamma_ray_test_composition, atomic_dataset, nuclide_name
):
    """
    Function to test if the mass-energy conservation is satisfied.
    Parameters
    ----------
    gamma_ray_test_composition: Function holding the composition.
    atomic_dataset: Tardis atomic-nuclear dataset
    nuclide_name: Name of the nuclide."""

    raw_isotopic_mass_fraction, cell_masses = gamma_ray_test_composition
    gamma_ray_lines = atomic_dataset.decay_radiation_data
    isotope_dict = create_isotope_dicts(raw_isotopic_mass_fraction, cell_masses)
    inventories_dict = create_inventories_dict(isotope_dict)
    total_decays = calculate_total_decays(inventories_dict, 1 * u.d)
    isotope_decay_df = create_isotope_decay_df(total_decays, gamma_ray_lines)

    grouped_isotope_df = isotope_decay_df.groupby(
        level=["shell_number", "isotope"]
    )

    parent_isotope_energy = (
        grouped_isotope_df.get_group((0, nuclide_name.replace("-", "")))[
            "energy_per_channel_keV"
        ].sum()
        * (u.keV).to(u.MeV)
        * u.MeV
    )

    neutrino_energy = 0.41 * u.MeV

    total_energy_actual = parent_isotope_energy + neutrino_energy

    c2 = const.c.to("cm/s") ** 2

    # calculate mass of 56Ni
    parent_isotope = rd.Nuclide(nuclide_name.replace("-", ""))
    parent_atomic_mass = parent_isotope.atomic_mass * (u.u).to(u.g) * u.g

    # calculate mass of 56Co
    daughter_isotope = parent_isotope.progeny()[0]

    daughter_atomic_mass = (
        rd.Nuclide(daughter_isotope).atomic_mass * (u.u).to(u.g) * u.g
    )

    Q = (parent_atomic_mass - daughter_atomic_mass) * c2 * u.erg.to(u.MeV)

    np.testing.assert_allclose(total_energy_actual.value, Q.value, rtol=0.01)


@pytest.mark.parametrize("nuclide_name", ["Ni-56", "Fe-52", "Cr-48"])
def test_activity(gamma_ray_test_composition, nuclide_name):
    """
    Function to test the decay of an atom in radioactivedecay with an analytical solution.
    Parameters
    ----------
    gamma_ray_test_composition: Function holding the composition.
    nuclide_name: Name of the nuclide.
    """
    # setup of decay test
    nuclide = rd.Nuclide(nuclide_name)
    t_half = nuclide.half_life() * u.s
    decay_constant = np.log(2) / t_half
    time_delta = 1.0 * u.s

    # calculating necessary values
    raw_isotopic_mass_fraction, cell_masses = gamma_ray_test_composition
    isotopic_masses = raw_isotopic_mass_fraction * cell_masses
    test_mass = isotopic_masses.loc[(nuclide.Z, nuclide.A), 0] * u.g
    isotope_dict = create_isotope_dicts(raw_isotopic_mass_fraction, cell_masses)
    inventories_dict = create_inventories_dict(isotope_dict)

    total_decays = calculate_total_decays(inventories_dict, time_delta)
    actual = total_decays.loc[
        (0, nuclide_name.replace("-", "")), "number_of_decays"
    ]

    isotope_mass = nuclide.atomic_mass * u.u
    number_of_atoms = (test_mass / isotope_mass).to(u.dimensionless_unscaled)
    expected = number_of_atoms * (1 - np.exp(-decay_constant * time_delta))

    npt.assert_allclose(actual, expected)


def test_total_energy_production(gamma_ray_test_composition, atomic_dataset):
    """
    Function to test the total energy production with equation 18 of Nadyozhin 1994. This is only for Ni56 now.
    Parameters
    ----------
    gamma_ray_test_composition: Function holding the composition.
    """
    time_start = 0.0 * u.d
    time_end = np.inf * u.d
    time_delta = (time_end - time_start).value

    gamma_ray_lines = atomic_dataset.decay_radiation_data
    raw_isotopic_mass_fraction, cell_masses = gamma_ray_test_composition
    isotope_dict = create_isotope_dicts(raw_isotopic_mass_fraction, cell_masses)
    inventories_dict = create_inventories_dict(isotope_dict)
    total_decays = calculate_total_decays(inventories_dict, time_delta)
    isotope_decay_df = create_isotope_decay_df(total_decays, gamma_ray_lines)
    taus, parents = get_taus(raw_isotopic_mass_fraction)

    ni56_nuclide = rd.Nuclide("Ni56")
    atomic_mass_unit = const.u.cgs.value
    tau_ni56 = taus["Ni-56"]
    tau_co56 = taus["Co-56"]

    ni56_mass_fraction = raw_isotopic_mass_fraction.loc[(28, 56)]

    ni_56_mass = sum(ni56_mass_fraction * cell_masses)
    ni_56_mass_solar = ni_56_mass / const.M_sun.cgs.value

    shell_number_0 = isotope_decay_df[
        isotope_decay_df.index.get_level_values("shell_number") == 0
    ]

    ni56 = shell_number_0[shell_number_0.index.get_level_values(1) == "Ni56"]
    ni_56_energy = ni56["energy_per_channel_keV"].sum()

    co_56 = shell_number_0[shell_number_0.index.get_level_values(1) == "Co56"]
    co_56_energy = co_56["energy_per_channel_keV"].sum()

    first_term = const.M_sun.cgs.value / (
        (tau_co56 - tau_ni56) * ni56_nuclide.atomic_mass * atomic_mass_unit
    )
    ni_term = (
        (ni_56_energy * (tau_co56 / tau_ni56 - 1) - co_56_energy)
        * first_term
        * KEV2ERG
    )

    co_term = co_56_energy * first_term * KEV2ERG

    expected = (
        ni_term
        * tau_ni56
        * (
            np.exp(-time_start.value / tau_ni56)
            - np.exp(-time_end.value / tau_ni56)
        )
        + co_term
        * tau_co56
        * (
            np.exp(-time_start.value / tau_co56)
            - np.exp(-time_end.value / tau_co56)
        )
    ) * ni_56_mass_solar

    ni56_df = isotope_decay_df[
        isotope_decay_df.index.get_level_values(1) == "Ni56"
    ]
    ni56_energy = ni56_df["decay_energy_erg"].sum()
    co_56_df = isotope_decay_df[
        isotope_decay_df.index.get_level_values(1) == "Co56"
    ]
    co56_energy = co_56_df["decay_energy_erg"].sum()
    actual = ni56_energy + co56_energy

    npt.assert_allclose(actual, expected)


def test_cumulative_decays(gamma_ray_test_composition, atomic_dataset):
    """
    Function to test that the total energy calculated from summing all the decays
    from the entire time range of simulation is the same as decay energy from individual
    time steps considering that after each time step the composition (mass fractions) changes.
    Tested for Ni56, Cr48, Fe52.
    Parameters
    ----------
    gamma_ray_simulation_state: Tardis simulation state
    atomic_dataset: Tardis atomic-nuclear dataset
    """

    time_start = 0.1 * u.d
    time_end = 100 * u.d
    time_steps = 3
    time_space = "linear"
    time_delta = (time_end - time_start).value

    gamma_ray_lines = atomic_dataset.decay_radiation_data
    raw_isotopic_mass_fraction, cell_masses = gamma_ray_test_composition
    isotope_dict = create_isotope_dicts(raw_isotopic_mass_fraction, cell_masses)
    inventories_dict = create_inventories_dict(isotope_dict)
    total_decays = calculate_total_decays(inventories_dict, time_delta)
    isotope_decay_df = create_isotope_decay_df(total_decays, gamma_ray_lines)

    times, effective_times = get_effective_time_array(
        time_start.value, time_end.value, time_space, time_steps
    )
    # total decay energy in the entire time range
    actual = isotope_decay_df["decay_energy_erg"].sum()

    # time evolve the decay energy
    evolve_decays_with_time = time_evolve_cumulative_decay(
        raw_isotopic_mass_fraction, cell_masses, gamma_ray_lines, times
    )
    expected = evolve_decays_with_time["decay_energy_erg"].sum()

    # This rtol is set since the decay energy is calculated with Fe52 (which has Mn-52m as a daughter)
    # The data is not available for Mn-52m in the decay_radiation_data
    # If we use any other isotope without a metastable state, the total decay energy matches exactly.
    npt.assert_allclose(actual, expected, rtol=1e-4)
