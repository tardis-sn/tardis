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



def test_datatable():
    """
    Creating datatable and checking various operations
    """
    df = pd.DataFrame( data = pd.np.random.randint(0, 100, (10, 5)) , columns = list('ABCED') )
    unit = pd.Series(['m', 's', 't', 'm/s', 'kg'])
    datatable = DataTable (df, units=['m', 's', 't', 'm/s', 'kg']) 
    
    """
    Assertions for units while copying, slicing and get unit function
    """
    assert unit.equals(datatable.units)
     
    datatable2 = datatable.copy()
    assert unit.equals(datatable2.units)
    
    columns = datatable.columns
    assert datatable[[columns[1], columns[2]]].units.equals(unit)
    
    assert datatable.get_unit("C") == 't'
    
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
    result = datatable.series_dot_product([1,2,3,4,5,6,7,8,9,10],"kg",datatable.columns)
    assert result.units.equals(pd.Series(['mkg','skg','tkg','m/skg','kgkg']))

    datatable.scalar_mul(2,"s",datatable.columns)
    assert datatable.units.equals(pd.Series(['ms','ss','ts','m/ss','kgs']))
    assert datatable.shape == datatable2.shape
