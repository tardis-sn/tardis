import numpy as np
import pandas as pd
from tardis.parser import parser
import pytest
import tempfile
import os

def test_parser_output():

    parser.parse_atomic_data("../si2_osc_kurucz")

#   check if the h5 file is formed
    assert os.path.exists(os.path.join(os.getcwd(),'atomic_data.h5'))

#   extract dataframes from h5 file
    df1 = pd.read_hdf('atomic_data.h5',key="energyLevels")
    df2 = pd.read_hdf('atomic_data.h5',key="oscillator_strengths")

#   delete the hdf file as it's no longer required
    os.remove('atomic_data.h5')

#   check the shapes of DataFrames
    assert df1.shape == (157,10)
    assert df2.shape == (4196,9)
