from typing import AsyncContextManager
from attr import has
import pytest
import numpy as np
import pandas as pd

from tardis.util import QDataFrame
from astropy import units as u



@pytest.fixture
def example_qdf():
    np.random.seed(1963)
    qdf = QDataFrame(index=pd.MultiIndex.from_arrays(np.random.choice([1, 2, 3, 4, 5], size=(10, 3)).T),
                    columns=['a', 'b', 'c'], 
                    data=np.random.random(size=(10, 3)), 
            units={'a':u.cm, 'b':u.s, 'c':u.joule})
    return qdf

@pytest.fixture
def example_df():
    noqdf = pd.DataFrame(index=pd.MultiIndex.from_arrays(np.random.choice([1, 2, 3, 4, 5], size=(10, 3)).T),
                 columns=['a', 'b', 'c'], 
                 data=np.random.random(size=(10, 3)))
    return noqdf
#qdf_from_noqdf = QDataFrame(noqdf, units={'a':u.cm, 'b':u.s, 'c':u.g})

def test_qdf_units(example_qdf):
    assert hasattr(example_qdf, 'units')
