import numpy as np
import pytest

pytestmark = pytest.mark.skip("Skipping tests due to old format")


def test_stimulated_emission_factor(stimulated_emission_factor):
    assert stimulated_emission_factor.shape == (253, 20)
    assert np.allclose(stimulated_emission_factor[250], 0.0010629785664215685)


def test_tau_sobolev(tau_sobolev):
    assert tau_sobolev.shape == (253, 20)
    assert np.allclose(tau_sobolev.iloc[565129], 3.0598663523617094e-05)


def test_beta_sobolev(beta_sobolev):
    assert beta_sobolev.shape == (253, 20)
    assert np.allclose(beta_sobolev[10][10], 1.671404577575537e-07)
