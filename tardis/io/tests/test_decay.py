import pandas as pd
import pytest
from astropy import units as u
from numpy.testing import assert_almost_equal

from tardis.model.matter.decay import IsotopicMassFraction


@pytest.fixture
def simple_abundance_model() -> IsotopicMassFraction:
    index = pd.MultiIndex.from_tuples(
        [(28, 56)], names=["atomic_number", "mass_number"]
    )
    return IsotopicMassFraction([[1.0, 1.0]], index=index)


def test_simple_decay(simple_abundance_model: IsotopicMassFraction) -> None:
    decayed_abundance = simple_abundance_model.calculate_decayed_mass_fractions(
        100 * u.day
    )
    assert_almost_equal(decayed_abundance.loc[26, 56][0], 0.55754786)
    assert_almost_equal(decayed_abundance.loc[26, 56][1], 0.55754786)
    assert_almost_equal(decayed_abundance.loc[27, 56][0], 0.44235126)
    assert_almost_equal(decayed_abundance.loc[27, 56][1], 0.44235126)
    assert_almost_equal(decayed_abundance.loc[28, 56][0], 1.10859709e-05)
    assert_almost_equal(decayed_abundance.loc[28, 56][1], 1.10859709e-05)


@pytest.fixture
def raw_abundance_simple() -> pd.DataFrame:
    abundances = pd.DataFrame([[0.2, 0.2], [0.1, 0.1]], index=[28, 30])
    abundances.index = abundances.index.rename("atomic_number")
    return abundances


def test_abundance_merge(
    simple_abundance_model: IsotopicMassFraction, raw_abundance_simple: pd.DataFrame
) -> None:
    decayed_df = simple_abundance_model.calculate_decayed_mass_fractions(100 * u.day)
    isotope_df = decayed_df.as_atoms()
    combined_df = decayed_df.merge(raw_abundance_simple, normalize=False)

    assert_almost_equal(
        combined_df.loc[28][0],
        raw_abundance_simple.loc[28][0] + isotope_df.loc[28][0],
    )
    assert_almost_equal(
        combined_df.loc[28][1],
        raw_abundance_simple.loc[28][1] + isotope_df.loc[28][1],
    )
    assert_almost_equal(combined_df.loc[30][1], raw_abundance_simple.loc[30][1])
    assert_almost_equal(combined_df.loc[26][0], isotope_df.loc[26][0])
