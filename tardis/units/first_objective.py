import pandas as pd
import numpy as np

class UDataFrame(pd.DataFrame):
  
    _metadata = ['units']

    @property
    def _constructor(self):
        def _c(*args, **kwargs):
            return UDataFrame(*args, **kwargs).__finalize__(self)
        return _c

    def __init__(self, *args, **kwargs):
        self.units = pd.Series(kwargs.pop("units", None))
        super().__init__(*args, **kwargs)

a = UDataFrame({"first_col": [1,2,3], "sec_col": [4, 5, 6]}, units=["cm", "mm"])

b = a.copy()

print(a, "Units for a:", a.units, " ",sep="\n"*2)

print("Units for b:", b.units,sep="\n"*2)
