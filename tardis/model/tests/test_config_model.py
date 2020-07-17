import numpy as np
import pandas as pd
import tardis
import os
from tardis.io.config_reader import Configuration
from tardis.model import Radial1DModel
import pytest

DATA_PATH = os.path.join(tardis.__path__[0],'model','tests','data')

@pytest.fixture(scope='module')
def config_hard_coded():
    configpath = os.path.join(DATA_PATH, 'config_hard_coded.yml')
    config = Configuration.from_yaml(configpath)
    config_hard_coded = Radial1DModel.from_config(config)
    return config_hard_coded 
 
def test_config_input_abundance_storage(config_hard_coded):

    #test abundance
    abundance_index = pd.Index([1,2,8,28],name = 'atomic_number')
    input_abundance = pd.DataFrame([[0.5,0.5,0.5,0.5,0.5], 
                                    [0.2,0.2,0.2,0.2,0.2],
                                    [0.25,0.25,0.25,0.25,0.25],
                                    [0.01,0.01,0.01,0.01,0.01]],
                                     index = abundance_index)

    pd.testing.assert_frame_equal(config_hard_coded.raw_abundance,input_abundance,check_names = True)
    
    #test isotopes
    arrays = [[28],[56]]
    isotope_index = pd.MultiIndex.from_arrays(
            arrays, names=["atomic_number", "mass_number"]
        )
    input_isotopes = pd.DataFrame([[0.04,0.04,0.04,0.04,0.04]],
                                  columns=np.arange(5),
                                  index = isotope_index)
    
    pd.testing.assert_frame_equal(config_hard_coded.raw_isotope_abundance,input_isotopes,check_names = True)
    

#test decay
def test_decay_dataframe(config_hard_coded):
    decay_index = pd.Index([1, 2, 8, 26, 27, 28],name = 'atomic_number')
    ref_abund = pd.DataFrame([
        [0.5,0.5,0.5,0.5,0.5], 
        [0.2,0.2,0.2,0.2,0.2],
        [0.25,0.25,0.25,0.25,0.25],
        [0.00096367034,0.00096367034,0.00096367034,0.00096367034,0.00096367034],
        [0.022980162,0.022980162,0.022980162,0.022980162,0.022980162],
        [0.026056167,0.026056167,0.026056167,0.026056167,0.026056167]
                             ], index = decay_index)

    pd.testing.assert_frame_equal(config_hard_coded.abundance,ref_abund,check_names = True, check_less_precise = True)
