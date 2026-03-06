import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from tardis.plasma.equilibrium.rates.radiative_rates import RadiativeRatesSolver


class DummyRadiationField:
    """A simple mock radiation field to isolate the RadiativeRatesSolver test."""
    def calculate_mean_intensity(self, nu):
        # Returns a dummy intensity array of 10.0 for all frequencies
        return np.full_like(nu, 10.0, dtype=float)

def test_radiative_rates_solver(regression_data):
    # 1. Setup specific MultiIndex required by the assertions
    index = pd.MultiIndex.from_tuples(
        [(2, 1, 0, 1), (2, 1, 1, 2)], # (atomic_number, ion_number, level_lower, level_upper)
        names=["atomic_number", "ion_number", "level_number_lower", "level_number_upper"]
    )

    # 2. Setup the DataFrame with required Einstein coefficients and frequency (nu)
    einstein_coefficients = pd.DataFrame(
        {
            "A_ul": [1.0e8, 2.0e8],
            "B_ul": [1.0e9, 2.0e9],
            "B_lu": [3.0e9, 4.0e9],
            "nu":   [1.0e15, 2.0e15]
        },
        index=index
    )
    # 3. Initialize solver
    solver = RadiativeRatesSolver(einstein_coefficients)

    # 4. Mock the radiation field and solve to get rates DataFrame
    radiation_field = DummyRadiationField()
    actual_rates_df = solver.solve(radiation_field)

    # 5. Save and compare with local regression data
    expected_rates_df = regression_data.sync_dataframe(actual_rates_df)
    assert_frame_equal(actual_rates_df, expected_rates_df)
