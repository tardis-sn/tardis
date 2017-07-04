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
