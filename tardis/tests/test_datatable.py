"""
=====================
TARDIS datatable test module
=====================
created on Mar 3, 2020
"""

import numpy as np
import pandas as pd
from datatable import DataTable
import pytest


@pytest.mark.parametrize ( " df , unit " ,
                         [pd.DataFrame( 
                                 data = pd.np.random.randint(0, 100, (10, 5)) , 
                                 columns = list('ABCED') ), 
                         pd.Series(['m', 's', 't', 'm/s', 'kg']) ] )


def test_datatable(df, unit):
    datatable = DataTable (df, units=['m', 's', 't', 'm/s', 'kg']) 

    assert unit.equals(datatable.units)

    datatable2 = datatable.copy()
    assert unit.equals(datatable2.units)

    columns = datatable.columns
    assert datatable[[columns[1], columns[2]]].units.equals(unit)

    assert datatable.get_unit("C") == 't'

    datatable[columns[0]] += 1
    assert unit.equals(datatable.units)

    datatable[columns[3]] *= 2
    assert unit.equals(datatable.units)
