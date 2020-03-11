import pandas as pd
import numpy as np
import UDataFrame.py

"""
Test of class UDataFrame by making an instance and checking of units preservation after copying and simple operations

"""

def test():

    units = pd.Series(["cm", "mm", "kg"])
    df = UDataFrame({"a": [1,2,3], "b": [4, 5, 6], "c": [7, 8, 9]}, units=units) 

    df2 = df.copy()
    assert units.equals(df2.units)

    df = np.log(df*2+1)
    assert units.equals(df.units)

test()
