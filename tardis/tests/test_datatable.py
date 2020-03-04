"""
=====================
TARDIS datatable test module
=====================
created on Mar 3, 2020
"""

import numpy as np
import pandas as pd
from tardis.datatable import DataTable
import pytest
import astropy.units as u


def test_datatable():
    """
    Creating datatable and checking various operations
    """
    df = pd.DataFrame( data = pd.np.random.randint(0, 100, (10, 5)) , columns = list('ABCED') )
    unit = pd.Series([u.meter, u.second, u.kg, u.meter/u.second, u.Ohm])
    datatable = DataTable (df, units=[u.meter, u.second, u.kg, u.meter/u.second, u.Ohm]) 
    
    """
    Assertions for units while copying, slicing and get unit function
    """
    assert unit.equals(datatable.units)
     
    datatable2 = datatable.copy()
    assert unit.equals(datatable2.units)
    
    columns = datatable.columns
    assert datatable[[columns[1], columns[2]]].units.equals(unit)
    
    assert datatable.get_unit("C") == u.kg
    
    """
    Assertion while arithmetic operations
    """
    datatable[columns[0]] += 1
    assert unit.equals(datatable.units)
    
    datatable[columns[3]] *= 2
    assert unit.equals(datatable.units)
    
    """
    Unit testing for scalar and series multiplication 
    """
    
    datatable.scalar_mul(2,u.second,datatable.columns)
    assert datatable.units.equals(pd.Series([u.meter*u.second, u.second*u.second, u.kg*u.second,
                                          u.meter*u.second/u.second, u.Ohm*u.second]))
    assert datatable.shape == datatable2.shape
    
    df = pd.DataFrame([[0, 1, -2, -1], [1, 1, 1, 1]])
    datatable3 = DataTable (df, units=[u.meter, u.second])
    s = pd.Series([1, 1, 2, 1])
    result = datatable3.dot(s, u.second, datatable3.columns)
    assert result.units.equals(pd.Series([u.meter*u.second, u.second*u.second]))
