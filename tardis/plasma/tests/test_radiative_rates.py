import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
from astropy import units as u
from tardis.plasma.equilibrium.rates.radiative_rates import RadiativeRatesSolver
from tardis.plasma.radiation_field.planck_rad_field import DilutePlanckianRadiationField

def test_radiative_rates_solver(regression_data):
    # 1. Setup specific MultiIndex required by the assertions
    index = pd.MultiIndex.from_tuples(
        [(2, 1, 0, 1), (2, 1, 1, 2)],
        names=["atomic_number", "ion_number", "level_number_lower", "level_number_upper"]
    )

    # 2. Setup the DataFrame with required Einstein coefficients
    einstein_coefficients = pd.DataFrame(
        {
            "A_ul": [1.0e8, 2.0e8],
            "B_ul": [1.0e9, 2.0e9],
            "B_lu": [3.0e9, 4.0e9],
            "nu":   [1.0e15, 2.0e15]
        },
        index=index
    )
    
    # 3. Setup a REAL Radiation Field using Astropy Units (NO MOCKING)
    t_rad = np.array([10000.0, 12000.0]) * u.K
    w = np.array([0.5, 0.4])
    real_radiation_field = DilutePlanckianRadiationField(temperature=t_rad, dilution_factor=w)

    # 4. Initialize solver and solve using the REAL radiation field
    solver = RadiativeRatesSolver(einstein_coefficients)
    actual_rates_df = solver.solve(real_radiation_field)

    # 5. Save and compare with local regression data
    expected_rates_df = regression_data.sync_dataframe(actual_rates_df)
    assert_frame_equal(actual_rates_df, expected_rates_df)
