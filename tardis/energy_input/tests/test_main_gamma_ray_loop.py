from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pandas as pd
from numpy.testing import assert_allclose

from tardis.energy_input.main_gamma_ray_loop import (
    calculate_electron_number_density,
)


def test_electron_density_includes_stable_elements() -> None:
    elemental_number_density = pd.DataFrame(
        {0: [2.0, 3.0, 4.0]}, index=pd.Index([6, 8, 28], name="atomic_number")
    )
    state = SimpleNamespace(
        v_inner=np.array([1.0]) * u.cm / u.s,
        v_outer=np.array([2.0]) * u.cm / u.s,
        volume=np.array([1.0]) * u.cm**3,
        calculate_elemental_number_density=lambda masses: elemental_number_density,
    )
    atom_data = SimpleNamespace(atom_data=SimpleNamespace(mass=pd.Series()))

    actual = calculate_electron_number_density(
        state,
        np.array([1.0]),
        np.array([1.0]),
        atom_data=atom_data,
    )

    expected = np.array(
        [(6 * 2 + 8 * 3 + 28 * 4) / (4 * np.pi / 3 * (2**3 - 1))]
    )
    assert_allclose(actual[:, 0], expected)
