import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.base import BasePlasma
from tardis.plasma.exceptions import IncompleteAtomicData
from tardis.plasma.properties.atomic import IonizationData

ionization_data_property = IonizationData(plasma_parent=BasePlasma)


def test_ionization_data_filter_atomic_property(atomic_dataset, regression_data):
    selected_atoms = [1, 2]
    actual_ionization_data = ionization_data_property._filter_atomic_property(
        atomic_dataset.ionization_data, selected_atoms
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
