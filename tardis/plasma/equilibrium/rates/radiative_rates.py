import numpy as np
import pandas as pd

from tardis.configuration.sorting_globals import SORTING_ALGORITHM


class RadiativeRatesSolver:
    """Calculate radiative bound-bound transition rates."""

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
        self.einstein_coefficients = einstein_coefficients.sort_index(
            kind=SORTING_ALGORITHM
        )

    def solve(
        self,
        radiation_field: object,
        beta_sobolevs: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
        """Calculate radiative transition rates.

        Parameters
        ----------
        radiation_field : RadiationField
            Radiation field used to calculate the line mean intensities.
        beta_sobolevs : pandas.DataFrame, optional
            Sobolev escape probabilities indexed by the Einstein transitions.
            When supplied, both radiative directions are multiplied by the
            corresponding escape probability.

        Returns
        -------
        pandas.DataFrame
            Transition rates indexed by source and destination states.
        """
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

        if beta_sobolevs is not None:
            beta_sobolevs = beta_sobolevs.reindex(
                self.einstein_coefficients.index
            )
            if beta_sobolevs.isna().any().any():
                raise ValueError(
                    "beta_sobolevs does not cover all Einstein transitions"
                )
            r_lu = r_lu.multiply(beta_sobolevs, axis=0)
            r_ul = r_ul.multiply(beta_sobolevs, axis=0)

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

        rates_df = rates_df.reset_index()

        # Add the new columns by duplicating the ion_number column
        rates_df["ion_number_source"] = rates_df["ion_number"]
        rates_df["ion_number_destination"] = rates_df["ion_number"]

        rates_df = rates_df.set_index(
            [
                "atomic_number",
                "ion_number",
                "ion_number_source",
                "ion_number_destination",
                "level_number_source",
                "level_number_destination",
            ]
        )

        return rates_df
