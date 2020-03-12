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
    datatable4 = datatable.copy()
    
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
    
    datatable.scalar_mul(2,u.second)
    assert datatable.units.equals(pd.Series([u.meter*u.second, u.second*u.second, u.kg*u.second,
                                          u.meter*u.second/u.second, u.Ohm*u.second]))
    assert datatable.shape == datatable2.shape
    
    df = pd.DataFrame([[0, 1, -2, -1], [1, 1, 1, 1]])
    datatable3 = DataTable (df, units=[u.meter, u.second])
    s = pd.Series([1, 1, 2, 1])
    result = datatable3.dot(s, u.second)
    assert result.units.equals(pd.Series([u.meter*u.second, u.second*u.second]))
    
    """
    Unit testing for add, replace and delete
    """
    datatable4.delete("E")
    units2 = pd.Series([u.meter, u.second, u.kg, u.Ohm])
    assert units2.equals(datatable4.units)
    
    units2 = pd.Series([u.meter, u.second, u.kg, u.Ohm,u.meter])
    datatable4.add(pd.Series([1,2,3,4,5,1,2,3,4,5]),u.meter,"F")
    assert datatable4.units.equals(units2)

    units2 = pd.Series([u.second, u.second, u.kg, u.Ohm,u.meter])
    datatable4.replace(pd.Series([1,2,3,4,5,1,2,3,4,5]),u.second,"A")
    assert datatable4.units.equals(units2)

    

test_datatable()
