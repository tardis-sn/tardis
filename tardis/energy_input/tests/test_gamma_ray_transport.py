import pytest
import numpy as np
import numpy.testing as npt
from tardis.model import SimulationState
from tardis.io.configuration import config_reader
from tardis.energy_input.gamma_ray_transport import calculate_shell_masses


# pytest.fixtures.
# Return the simulation object.


@pytest.fixture()
def config():
    return config_reader.Configuration.from_yaml(
        "/Users/anirbandutta/Projects/gamma_ray_tardis/tardis_configv1_density_exponential_nebular.yml"
    )


def test_calculate_shell_masses(config):
    import astropy.units as u

    config.model.structure.velocity.start = 1.0 * u.km / u.s
    config.model.structure.density.rho_0 = 5.0e2 * u.g / (u.cm**3)
    model = SimulationState.from_config(config)
    volume = 4.2006589e21 * (u.cm**3)
    density = 3.3848916e9 * u.g / (u.cm**3)

    shell_masses = calculate_shell_masses(model)
    print(shell_masses)
    actual = shell_masses[0].value
    desired = (volume * density).value
    npt.assert_almost_equal(actual, desired)
