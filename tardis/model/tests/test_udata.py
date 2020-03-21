import os
import pandas as pd
import pytest
import astropy.units as u
from ..udata import UDataFrame


@pytest.fixture(scope="module")
def units_series():
    return pd.Series([u.meter, None, u.kg, None, u.Ohm])

@pytest.fixture(scope="module")
def units_list():
    return [u.meter, None, u.kg, None, u.Ohm]

@pytest.fixture(scope="module")
def sample_data():
    return pd.DataFrame(data = pd.np.random.randint(0, 50, (10, 5)),
                        columns=['A','B','C','D','E'])

@pytest.fixture(scope="module")
def sample_udata(sample_data, units_list):
    return UDataFrame(sample_data, units=units_list)

def test_creation_utable(sample_udata, units_list, units_series):
    for i in range(len(units_series)):
        assert sample_udata.units[i] == units_series[i]
        assert sample_udata.units[i] == units_list[i]

    assert sample_udata.get_unit('A') == u.meter
    with pytest.raises(ValueError):
        sample_udata.get_unit('F')
    
    sample_udata.set_unit('B', u.second)
    assert sample_udata.get_unit('B') == u.second

def test_copy(sample_udata):
    new_data = sample_udata.copy()
    assert new_data.units.equals(sample_udata.units)
