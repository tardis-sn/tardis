import os

import numpy as np
import pandas as pd
import pytest
from astropy import units as u

import tardis
from tardis.plasma.properties import *
from tardis.io.atom_data import AtomData

# INPUTS

@pytest.fixture
def atomic_data(regression_data):
    atom_data_path = (
    regression_data.regression_data_path
        / "atom_data"
        / "new_kurucz_cd23_chianti_H_He.h5" 
    )
    return AtomData.from_hdf(atom_data_path)

