import numpy as np

def test_stimulated_emission_factor(stimulated_emission_factor):
    assert stimulated_emission_factor.shape == (247,20)
    assert np.allclose(stimulated_emission_factor[244], 0.0010629785664215685)

def test_tau_sobolev(tau_sobolev):
    assert tau_sobolev.shape == (247,20)
    assert np.allclose(tau_sobolev.ix[533329], 3.04484722689e-05)

def test_beta_sobolev(beta_sobolev):
    assert beta_sobolev.shape == (247,20)
    assert np.allclose(beta_sobolev.iloc[10,10], 1.671404577575537e-07)
