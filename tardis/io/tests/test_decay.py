import pytest
import pandas as pd

from tardis.io.decay import IsotopeAbundances

@pytest.fixture
def simple_abundance_model():
    index = pd.MultiIndex.from_tuples([(28, 56)],
                                      names=['atomic_number', 'mass_number'])
    return IsotopeAbundances([[1.0, 1.0]], index=index)


def test_simple_decay(simple_abundance_model):
    1/0
