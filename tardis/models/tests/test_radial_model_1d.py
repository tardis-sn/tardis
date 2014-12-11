import pytest

import numpy as np

import numpy.testing as nptesting

from astropy import units as u

from tardis.models.model1d import Radial1D, HomologousRadial1D


@pytest.fixture
def radius():
    return np.arange(1000, 11000, 1000) * u.km

@pytest.fixture
def time():
    return 10 * u.day

@pytest.fixture
def time0():
    return 2 * u.s

@pytest.fixture
def density0():
    return np.linspace(20, 1, 20) * u.g / u.cm**3

@pytest.fixture
def velocity():
    return np.linspace(10000, 30000, 20) * u.km / u.s

@pytest.fixture
def simple_rh_model(velocity, time, time0, density0):
    return HomologousRadial1D(velocity, time, time0, density0)


def test_simple_radial_model1d_1(radius):

    model = Radial1D(radius)

    assert model.radius.unit == u.cm


def test_simple_homologous_model(velocity, time, time0):
    model = HomologousRadial1D(velocity, time, time0)

    assert model.radius.unit == u.cm
    assert model.velocity.unit == u.cm / u.s

    nptesting.assert_allclose(model.volume[0],
                              4. / 3. * np.pi * model.r_inner[0]**3)


def test_len_check():
    with pytest.raises(AttributeError):
        model = HomologousRadial1D(np.arange(1, 10) * u.km/u.s, time=10*u.day,
                               density0=np.arange(1, 20), time0=5)


def test_homologous_density(simple_rh_model):
    model = simple_rh_model

    nptesting.assert_allclose(model.density, (model.time0/model.time)**3
                              * model.density0)
    nptesting.assert_allclose(model.mass, model.density * model.volume)