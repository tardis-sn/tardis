import pytest

import numpy as np

from astropy import units as u

from tardis.models.model1d import RadialModel1D


def test_simple_radial_model1d_1():
    radius = np.arange(1000, 11000, 1000) * u.km
    density = np.linspace(1e8, 1e9, 10) * u.Msun / u.AU**3
    test_model = RadialModel1D(radius, density, None, None)

    assert test_model.radius.unit == u.cm
    assert test_model.density.unit == u.g / u.cm**3
