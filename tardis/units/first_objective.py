"""
==========================================

UDataFrame units extantion of pd.DataFrame 

==========================================
"""

import pandas as pd
import numpy as np

class UDataFrame(pd.DataFrame):
    """
    Single subclassing of pd.DataFrame.
    This class variable tells Pandas the name of the attributes
    that are to be ported over to derivative DataFrames.
    """
    _metadata = ['units']

    @property
    def _constructor(self):
        '''
        This is the key to letting Pandas know how to keep
        derivative `UDataFrame` the same type as my. It should
        be enough to return the name of the Class.
        There is a method named `__finalize__` that grabs these attributes
        and assigns them to newly created `UDataFrame`
        Workaround for https://github.com/pandas-dev/pandas/issues/13208
        '''
        def _c(*args, **kwargs):
            return UDataFrame(*args, **kwargs).__finalize__(self)
        return _c

    def __init__(self, *args, **kwargs):
        """
        Grab the keyword argument that is supposed to be units
        """
        self.units = pd.Series(kwargs.pop("units", None))
        super().__init__(*args, **kwargs)
