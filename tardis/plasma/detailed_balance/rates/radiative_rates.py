import numpy as np
import pandas as pd


class RadiativeRatesSolver:
    einstein_coefficients: pd.DataFrame

    def __init__(self, einstein_coefficients):
        # Ensuring the right columns are present
        assert einstein_coefficients.index.names == [
            "atomic_number",
            "ion_number",
            "level_number_lower",
            "level_number_upper",
        ]
        assert {"A_ul", "B_ul", "B_lu", "nu"} - set(
            einstein_coefficients.columns
        ) == set()

        assert np.all(
            einstein_coefficients.index.get_level_values("level_number_lower")
            < einstein_coefficients.index.get_level_values("level_number_upper")
        )
        self.einstein_coefficients = einstein_coefficients.sort_index()

    def solve(self, radiation_field):
        mean_intensity = radiation_field.calculate_mean_intensity(
            self.einstein_coefficients.nu.values
        )
        mean_intensity_df = pd.DataFrame(
            data=mean_intensity, index=self.einstein_coefficients.index
        )

        # r_lu = B_lu * J_nu
        r_lu = mean_intensity_df.multiply(
            self.einstein_coefficients.B_lu, axis=0
        )

        # r_ul = B_ul * J_nu + A_ul
        r_ul = mean_intensity_df.multiply(
            self.einstein_coefficients["B_ul"], axis=0
        )
        r_ul = r_ul.add(self.einstein_coefficients["A_ul"], axis=0)

        # swapping as source is upper and destination is lower
        r_ul.index = r_ul.index.swaplevel(
            "level_number_lower", "level_number_upper"
        )

        rates_df = pd.concat([r_lu, r_ul])
        rates_df.index.names = [
            "atomic_number",
            "ion_number",
            "level_number_source",
            "level_number_destination",
        ]
        return rates_df
