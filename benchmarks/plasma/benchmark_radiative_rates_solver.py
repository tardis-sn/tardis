import numpy as np
import pandas as pd

from tardis.plasma.equilibrium.rates.radiative_rates import (
    RadiativeRatesSolver,
)


class DummyRadiationField:
    def calculate_mean_intensity(self, nu):
        return np.ones_like(nu)


def make_einstein_coeff_df(n_lines):
    index = pd.MultiIndex.from_tuples(
        [
            (1, 0, i, i + 1)
            for i in range(n_lines)
        ],
        names=[
            "atomic_number",
            "ion_number",
            "level_number_lower",
            "level_number_upper",
        ],
    )

    return pd.DataFrame(
        {
            "A_ul": np.ones(n_lines),
            "B_ul": np.ones(n_lines),
            "B_lu": np.ones(n_lines),
            "nu": np.linspace(1e14, 5e14, n_lines),
        },
        index=index,
    )


class BenchmarkRadiativeRatesSolver:
    """
    ASV benchmark for RadiativeRatesSolver.solve()

    Measures runtime performance as number of transitions increases.
    """

    params = [50, 100, 200]
    param_names = ["n_lines"]

    def setup(self, n_lines):
        self.df = make_einstein_coeff_df(n_lines)
        self.solver = RadiativeRatesSolver(self.df)
        self.rad_field = DummyRadiationField()

    def time_solve(self, n_lines):
        self.solver.solve(self.rad_field)