import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis import constants as const
from tardis.plasma.base import BasePlasma
from tardis.plasma.exceptions import IncompleteAtomicData
from tardis.plasma.properties.atomic import IonizationData, Levels

ionization_data_property = IonizationData(plasma_parent=BasePlasma)
levels_property = Levels(plasma_parent=BasePlasma)


def test_beta_rad(t_rad, beta_rad):
    np.testing.assert_allclose(
        beta_rad,
        1 / (const.k_B.cgs.value * t_rad),
        atol=0,
        rtol=1e-15,
    )


def test_ionization_data_calculate_atomic_property(atomic_dataset, regression_data):
    selected_atoms = [1, 2]
    actual_ionization_data = ionization_data_property.calculate(
        atomic_dataset, selected_atoms
    )
    expected_ionization_data = regression_data.sync_dataframe(
        actual_ionization_data, key="ionization_data"
    )
    pdt.assert_series_equal(
        actual_ionization_data, expected_ionization_data, atol=0, rtol=1e-15
    )


@pytest.mark.xfail(raises=IncompleteAtomicData)
def test_ionization_data_incomplete_atomic_data():
    selected_atoms = [1, 2]
    index = pd.MultiIndex.from_tuples(
        [(1, 1), (2, 1)],
        names=["atomic_number", "ion_number"],
    )
    ionization_data = pd.Series(
        [1.0, 2.0], index=index, name="ionization_energy"
    )

    ionization_data_property._filter_atomic_property(
        ionization_data, selected_atoms
    )

def test_levels_calculate(atomic_dataset, regression_data):
    selected_atoms = [1, 2]
    levels_index, _, _, _ = levels_property.calculate(
        atomic_dataset, selected_atoms
    )
    actual_levels_calculate = atomic_dataset.levels.loc[levels_index]
    expected_levels_calculate = regression_data.sync_dataframe(
        actual_levels_calculate, key="filtered_levels"
    )
    pdt.assert_frame_equal(
        actual_levels_calculate, expected_levels_calculate, atol=0, rtol=1e-15
    )
