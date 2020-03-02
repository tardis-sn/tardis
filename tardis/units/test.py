import pandas as pd
import numpy as np
import UDataFrame

"""
Test of class UDataFrame by making an instance and checking of units preservation after copying
"""

a = UDataFrame({"first_col": [1,2,3], "sec_col": [4, 5, 6]}, units=["cm", "mm"])

b = a.copy()

print(a, "Units for a:", a.units, " ",sep="\n"*2)

print("Units for b:", b.units,sep="\n"*2)
