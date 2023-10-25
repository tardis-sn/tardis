import pytest
import numpy as np
import numpy.testing as npt

from tardis.model import SimulationState
from tardis.io.configuration import config_reader
from tardis.energy_input.gamma_ray_transport import (
    calculate_shell_masses,
    create_isotope_dicts,
    get_all_isotopes,
    create_inventories_dict,
    calculate_total_decays,
)
import astropy.units as u
import astropy.constants as c

# pytest.fixtures.
# Return the simulation object.


@pytest.fixture()
def config():
    return config_reader.Configuration.from_yaml(
        "/Users/anirbandutta/Projects/gamma_ray_tardis/tardis_configv1_density_exponential_nebular.yml"
    )


@pytest.fixture(scope="module", autouse=True)
def simulation_setup(config):
    config.model.structure.velocity.start = 1.0 * u.km / u.s
    config.model.structure.density.rho_0 = 5.0e2 * u.g / (u.cm**3)
    config.supernova.time_explosion = 1.0 * (u.d)
    model = SimulationState.from_config(config)
    return model


def test_activity(simulation_setup):
    """
    ni_mass : mass of 56Ni in grams.
    """
    model = simulation_setup
    t_half = 6.075 * u.d.to(u.s)  # day
    decay_constant = 0.693 / t_half  # sec^-1
    time_delta = 80.0 * u.d.to(u.s)  # day
    shell_masses = calculate_shell_masses(model)
    raw_isotope_abundance = model.raw_isotope_abundance
    raw_isotope_abundance_mass = raw_isotope_abundance.apply(
        lambda x: x * shell_masses, axis=1
    )
    ni_mass = raw_isotope_abundance_mass.loc[28, 56][0]
    iso_dict = create_isotope_dicts(raw_isotope_abundance, shell_masses)
    inv_dict = create_inventories_dict(iso_dict)
    total_decays = calculate_total_decays(inv_dict, time_delta)
    simulated_activity = total_decays[0][28, 56]["Ni-56"]

    isotopic_mass = 55.942132 * (u.g)  # g
    number_of_moles = ni_mass * (u.g) / isotopic_mass
    number_of_atoms = number_of_moles * c.N_A

    # Use the radioactivity formula
    # N = N_0 e-(lambda * t)
    N = number_of_atoms.value * np.exp(-decay_constant * time_delta)

    actual_activity = decay_constant * N
    npt.assert_almost_equal(actual_activity, simulated_activity)


@pytest.mark.xfail(reason="To be implemented")
def test_calculate_shell_masses(simulation_setup):
    model = simulation_setup
    volume = 4.2006589e21 * (u.cm**3)
    density = 3.3848916e9 * u.g / (u.cm**3)

    shell_masses = calculate_shell_masses()
    actual = shell_masses[0].value
    desired = (volume * density).value
    npt.assert_almost_equal(actual, desired)
