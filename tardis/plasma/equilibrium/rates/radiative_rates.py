import numpy as np
import pandas as pd

from tardis.configuration.sorting_globals import SORTING_ALGORITHM


class RadiativeRatesSolver:
    """Calculate bound-bound Einstein transition rates."""

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
        mean_intensity: pd.DataFrame,
        beta_sobolev: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
        """Calculate line rates from fixed intensities and escape probabilities.

        Parameters
        ----------
        mean_intensity : pandas.DataFrame
            Mean intensity blueward of each line. Rows follow the Einstein-line
            index and columns are shells.
        beta_sobolev : pandas.DataFrame, optional
            Candidate Sobolev escape probabilities aligned with
            ``mean_intensity``. When omitted, all escape probabilities are one.

        Returns
        -------
        pandas.DataFrame
            Upward and downward rates in source--destination index convention.
        """
        mean_intensity_df = mean_intensity.loc[self.einstein_coefficients.index]

        if beta_sobolev is None:
            beta_sobolev_array = np.ones(mean_intensity_df.shape)
        else:
            beta_sobolev_array = beta_sobolev.loc[
                self.einstein_coefficients.index, mean_intensity_df.columns
            ].to_numpy()

        if (
            not np.isfinite(mean_intensity_df.to_numpy()).all()
            or not np.isfinite(beta_sobolev_array).all()
        ):
            raise ValueError(
                "Mean intensities and Sobolev escape probabilities must be finite."
            )

        # r_lu = B_lu * J_nu
        r_lu = mean_intensity_df.multiply(
            self.einstein_coefficients.B_lu, axis=0
        )
        r_lu *= beta_sobolev_array

        # r_ul = B_ul * J_nu + A_ul
        r_ul = mean_intensity_df.multiply(
            self.einstein_coefficients["B_ul"], axis=0
        )
        r_ul = r_ul.add(self.einstein_coefficients["A_ul"], axis=0)
        r_ul *= beta_sobolev_array

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
        rates_df.columns = mean_intensity_df.columns.copy()

        return rates_df
