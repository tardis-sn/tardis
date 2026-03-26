import astropy.units as u
import numpy as np
import pandas.testing as pdt
import pytest

from tardis.plasma.equilibrium.rates.radiative_rates import RadiativeRatesSolver
from tardis.plasma.radiation_field.planck_rad_field import PlanckianRadiationField
from tardis.io.atom_data import AtomData

@pytest.fixture(scope="class")
def einstein_coefficients_df():
    atom_data = AtomData.from_hdf('kurucz_cd23_chianti_H_He_latest.h5')
    lines_df = atom_data.lines
    radiative_transitions = lines_df.loc[(1,0, slice(None), slice(None)), :]
    return radiative_transitions


@pytest.fixture(scope="class")
def mock_radiation_field():
    """Fixture for mock radiation field."""
    temperature = [10000, 20000]* u.K
    return PlanckianRadiationField(temperature=temperature)

def test_radiative_rate_solver_init(einstein_coefficients_df,regression_data):
    solver = RadiativeRatesSolver(einstein_coefficients_df)
    actual_einstein_coeffs = solver.einstein_coefficients
    expected_einstein_coeffs = regression_data.sync_dataframe(
        actual_einstein_coeffs, key="einstein_coeffs")
    pdt.assert_frame_equal(actual_einstein_coeffs,expected_einstein_coeffs)

def test_radiative_rate_solver_solve(einstein_coefficients_df, mock_radiation_field, regression_data):
    solver = RadiativeRatesSolver(einstein_coefficients_df)
    actual_radiative_rates = solver.solve(mock_radiation_field)
    expected_radiative_rates = regression_data.sync_dataframe(
        actual_radiative_rates, key="radiative_rates"
    )
    pdt.assert_frame_equal(actual_radiative_rates,expected_radiative_rates)
