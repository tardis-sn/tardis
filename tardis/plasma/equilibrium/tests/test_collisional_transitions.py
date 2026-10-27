from pathlib import Path

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy import units as u
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis.io.atom_data import AtomData
from tardis.plasma.equilibrium.rates import ThermalCollisionalRateSolver


@pytest.fixture
def cmfgen_atomic_data(tardis_regression_path: Path) -> AtomData:
    atom_data = AtomData.from_hdf(
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )
    atom_data.prepare_atom_data([1], "macroatom", [(1, 0)], [(1, 0)])
    atom_data.yg_data.columns = list(atom_data.collision_data_temperatures)
    return atom_data


def test_iip_cmfgen_collisional_strengths(
    cmfgen_atomic_data: AtomData,
    regression_data: RegressionData,
) -> None:
    """Verify IIP collision-strength interpolation at tabulated temperatures."""
    actual = ThermalCollisionalRateSolver(
        cmfgen_atomic_data.levels,
        cmfgen_atomic_data.lines,
        cmfgen_atomic_data.collision_data_temperatures,
        cmfgen_atomic_data.yg_data,
        collision_strengths_type="cmfgen",
    ).calculate_collision_strengths(
        cmfgen_atomic_data.collision_data_temperatures * u.K
    )
    expected = regression_data.sync_dataframe(
        pd.DataFrame(actual.to_numpy()), key="allclose_0"
    )
    npt.assert_allclose(actual.to_numpy(), expected.to_numpy(), rtol=1e-8, atol=0.0)


def test_thermal_collision_rates_against_iip(
    cmfgen_atomic_data: AtomData,
    regression_data: RegressionData,
) -> None:
    """Verify equilibrium collisional rates against IIP plasma properties."""
    transition_index = (1, 0, slice(None), slice(None))
    lines = cmfgen_atomic_data.lines.loc[transition_index, :]
    strengths = cmfgen_atomic_data.yg_data.loc[transition_index, :]
    actual = ThermalCollisionalRateSolver(
        cmfgen_atomic_data.levels,
        lines,
        cmfgen_atomic_data.collision_data_temperatures,
        strengths,
        collision_strengths_type="cmfgen",
        collisional_strength_approximation="regemorter",
    ).solve(np.array([10000.0, 20000.0]) * u.K)
    expected = regression_data.sync_dataframe(actual, key="frame_0")
    pdt.assert_frame_equal(
        actual,
        expected,
        check_names=False,
        check_column_type=False,
        atol=0.0,
        rtol=2e-5,
    )
