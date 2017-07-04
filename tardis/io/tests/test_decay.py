import pytest
import pandas as pd

from tardis.io.decay import IsotopeAbundances
from numpy.testing import assert_almost_equal

@pytest.fixture
def simple_abundance_model():
    index = pd.MultiIndex.from_tuples([(28, 56)],
                                      names=['atomic_number', 'mass_number'])
    return IsotopeAbundances([[1.0, 1.0]], index=index)

def test_simple_decay(simple_abundance_model):
    decayed_abundance = simple_abundance_model.decay(100)
    assert_almost_equal(decayed_abundance.ix[26, 56][0], 0.55752)
    assert_almost_equal(decayed_abundance.ix[26, 56][1], 0.55752)
    assert_almost_equal(decayed_abundance.ix[27, 56][0], 0.4423791)
    assert_almost_equal(decayed_abundance.ix[27, 56][1], 0.4423791)
    assert_almost_equal(decayed_abundance.ix[28, 56][0], 1.1086e-05)
    assert_almost_equal(decayed_abundance.ix[28, 56][1], 1.1086e-05)

@pytest.fixture
def raw_abundance_simple():
    abundances = pd.DataFrame([[0.2, 0.2], [0.1, 0.1]], index=[28, 30])
    abundances.index.rename('atomic_number', inplace=True)
    return abundances

def test_abundance_merge(simple_abundance_model, raw_abundance_simple):
    decayed_df = simple_abundance_model.decay(100)
    merged_df = simple_abundance_model.as_atomic_numbers(raw_abundance_simple, 100, normalize=False)

    assert_almost_equal(merged_df.loc[28][0], raw_abundance_simple.loc[28][0] + decayed_df.loc[28, 56][0])
    assert_almost_equal(merged_df.loc[28][1], raw_abundance_simple.loc[28][1] + decayed_df.loc[28, 56][1])
    assert_almost_equal(merged_df.loc[30][1], raw_abundance_simple.loc[30][1])
    assert_almost_equal(merged_df.loc[26][0], decayed_df.loc[26, 56][0])