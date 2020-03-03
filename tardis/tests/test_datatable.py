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


@pytest.mark.parametrize ( (" df" , "unit ") ,
                         [pd.DataFrame( 
                                 data = pd.np.random.randint(0, 100, (10, 5)) , 
                                 columns = list('ABCED') ), 
                         pd.Series(['m', 's', 't', 'm/s', 'kg']) ] )

def test_datatable(df, unit):
    assert 1==1
    
