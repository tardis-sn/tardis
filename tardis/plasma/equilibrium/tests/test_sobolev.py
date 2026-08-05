import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from tardis.opacities.tau_sobolev import (
    calculate_beta_sobolev,
    calculate_sobolev_line_opacity,
)
from tardis.plasma.equilibrium.sobolev import CoupledSobolevPopulationSolver
from tardis.plasma.properties.radiative_properties import (
    calculate_stimulated_emission_factor,
)

LEVEL_INDEX = pd.MultiIndex.from_tuples(
    [(1, 0, 0), (1, 0, 1)],
    names=["atomic_number", "ion_number", "level_number"],
)
LINE_INDEX = pd.MultiIndex.from_tuples(
    [(1, 0, 0, 1)],
    names=[
        "atomic_number",
        "ion_number",
        "level_number_lower",
        "level_number_upper",
    ],
)
LINES = pd.DataFrame(
    {"wavelength_cm": [1e-5], "f_lu": [1000.0]}, index=LINE_INDEX
)
G = pd.Series([1.0, 1.0], index=LEVEL_INDEX)
METASTABILITY = pd.Series([False, False], index=LEVEL_INDEX)


def solve_population(beta_sobolev: pd.DataFrame) -> pd.DataFrame:
    """Return a deterministic population whose excitation responds to beta."""
    upper_population = 0.2 + 0.4 * beta_sobolev.iloc[0, 0]
    return pd.DataFrame(
        [[1.0], [upper_population]], index=LEVEL_INDEX, columns=[0]
    )


def make_solver(
    time_explosion: u.Quantity = 1 * u.day,
) -> CoupledSobolevPopulationSolver:
    return CoupledSobolevPopulationSolver(
        solve_population,
        LINES,
        time_explosion,
        G,
        [0],
        [1],
        METASTABILITY,
        tolerance=1e-7,
    )


@pytest.mark.parametrize("upper_population", [0.1, 0.8])
def test_sobolev_fixed_point_converges_from_positive_states(
    upper_population: float,
) -> None:
    initial_level_number_density = pd.DataFrame(
        [[1.0], [upper_population]], index=LEVEL_INDEX, columns=[0]
    )

    result = make_solver().solve(initial_level_number_density)

    level_number_density, _, tau_sobolevs, beta_sobolev = result
    assert np.all(level_number_density.to_numpy() > 0.0)
    np.testing.assert_allclose(
        beta_sobolev.to_numpy(),
        calculate_beta_sobolev(tau_sobolevs).to_numpy(),
    )


