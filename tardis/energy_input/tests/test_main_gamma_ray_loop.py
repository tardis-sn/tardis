from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis.energy_input.main_gamma_ray_loop import (
    calculate_scaled_elemental_number_density,
)


@pytest.mark.parametrize("atomic_number_exponent", [1, 5])
def test_elemental_number_density_includes_stable_elements(
    atomic_number_exponent: int,
) -> None:
    elemental_number_density = pd.DataFrame(
        {0: [2.0, 3.0, 4.0]}, index=pd.Index([6, 8, 28], name="atomic_number")
    )
    state = SimpleNamespace(
        v_inner=np.array([1.0]) * u.cm / u.s,
        v_outer=np.array([2.0]) * u.cm / u.s,
        volume=np.array([1.0]) * u.cm**3,
        calculate_elemental_number_density=lambda masses: (
            elemental_number_density
        ),
    )
    atom_data = SimpleNamespace(atom_data=SimpleNamespace(mass=pd.Series()))

    actual = calculate_scaled_elemental_number_density(
        state,
        np.array([1.0]),
        np.array([1.0]),
        atom_data=atom_data,
        atomic_number_exponent=atomic_number_exponent,
    )

    expected = np.array(
        [
            (
                6**atomic_number_exponent * 2
                + 8**atomic_number_exponent * 3
                + 28**atomic_number_exponent * 4
            )
            / (4 * np.pi / 3 * (2**3 - 1))
        ]
    )
    assert_allclose(actual[:, 0], expected)


def test_scaled_elemental_number_density_regression(
    regression_data: RegressionData,
) -> None:
    elemental_number_density = pd.DataFrame(
        {0: [2.0, 3.0, 4.0]}, index=pd.Index([6, 8, 28], name="atomic_number")
    )
    state = SimpleNamespace(
        v_inner=np.array([1.0]) * u.cm / u.s,
        v_outer=np.array([2.0]) * u.cm / u.s,
        volume=np.array([1.0]) * u.cm**3,
        calculate_elemental_number_density=lambda masses: (
            elemental_number_density
        ),
    )
    atom_data = SimpleNamespace(atom_data=SimpleNamespace(mass=pd.Series()))
    actual = np.column_stack(
        [
            calculate_scaled_elemental_number_density(
                state,
                np.array([1.0]),
                np.array([1.0]),
                atom_data=atom_data,
                atomic_number_exponent=atomic_number_exponent,
            )[:, 0]
            for atomic_number_exponent in (1, 5)
        ]
    )

    assert_allclose(actual, regression_data.sync_ndarray(actual))
