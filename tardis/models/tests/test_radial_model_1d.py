import pytest

import numpy as np

import numpy.testing as nptesting

from astropy import units as u

from tardis.models.model1d import Radial1D, RadialHomologous1D


@pytest.fixture
def radius():
    return np.arange(1000, 11000, 1000) * u.km

@pytest.fixture
def time():
    return 10 * u.day

@pytest.fixture
def velocity():
    return np.linspace(10000, 30000, 20) * u.km / u.s
def test_simple_radial_model1d_1(radius):

    model = Radial1D(radius)

    assert model.radius.unit == u.cm
    #assert model.density.unit == u.g / u.cm**3


#    nptesting.assert_allclose(test_model.mass[0],  2.488679615529453e+26 * u.g)


def test_simple_homologous_model(velocity, time):
    model = RadialHomologous1D(velocity, time)

    assert model.radius.unit == u.cm
    assert model.velocity.unit == u.cm / u.s

    nptesting.assert_allclose(model.volume[0],
                              4. / 3. * np.pi * model.r_inner[0]**3)


