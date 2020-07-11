import numpy as np
import pandas as pd
import numpy.testing as npt
import tardis
import os
from astropy import units as u
from tardis.io.config_reader import Configuration
from tardis.model import Radial1DModel
import pytest

DATA_PATH = os.path.join(tardis.__path__[0],'model','tests','data')

@pytest.fixture(scope="module", params=[('config_csvy_full.yml','config_csvy_full_old.yml'),
                                        ('config_csvy_nocsv_branch85.yml','config_csvy_nocsv_branch85_old.yml'),
                                        ('config_csvy_nocsv_uniform.yml','config_csvy_nocsv_uniform_old.yml'),
                                        ('config_csvy_nocsv_powerlaw.yml','config_csvy_nocsv_powerlaw_old.yml'),
                                        ('config_csvy_nocsv_exponential.yml','config_csvy_nocsv_exponential_old.yml'),
                                        ('config_csvy_full_rad.yml','config_csvy_full_rad_old.yml')])
def filename(request):
    fn_tup = request.param
    fn = os.path.join(DATA_PATH, fn_tup[0])
    fnold = os.path.join(DATA_PATH, fn_tup[1])
    return fn, fnold

def test_compare_models(filename):
    fn, fnold = filename
    tardis_config = Configuration.from_yaml(fn)
    tardis_config_old = Configuration.from_yaml(fnold)
    csvy_model = Radial1DModel.from_csvy(tardis_config)
    config_model = Radial1DModel.from_config(tardis_config_old)
    csvy_model_props = csvy_model.get_properties().keys()
    config_model_props = config_model.get_properties().keys()
    npt.assert_array_equal(csvy_model_props, config_model_props)
    for prop in config_model_props:
        csvy_model_val = csvy_model.get_properties()[prop]
        config_model_val = config_model.get_properties()[prop]
        if prop == 'homologous_density':
            npt.assert_array_almost_equal(csvy_model_val.density_0.value, config_model_val.density_0.value)
            npt.assert_array_almost_equal(csvy_model_val.time_0.value, config_model_val.time_0.value)
        else:
            if hasattr(config_model_val, 'value'):
                config_model_val = config_model_val.value
                csvy_model_val = csvy_model_val.value
            npt.assert_array_almost_equal(csvy_model_val, config_model_val)

@pytest.fixture(scope='module')
def csvy():
    csvypath = os.path.join(DATA_PATH, 'config_v_filter.yml')
    config = Configuration.from_yaml(csvypath)
    csvy = Radial1DModel.from_csvy(config)
    return csvy

def test_csvy_input_abundance_storage(csvy):

    #test abundance
    abundance_index = pd.Index([1,2],name = 'atomic_number')
    input_abundance = pd.DataFrame([[0.0,0.4],[0.99,0.58]],index = abundance_index)
    pd.testing.assert_frame_equal(csvy.raw_abundance,input_abundance,check_names = True)
    
    #test isotopes
    arrays = [[28],[56]]
    isotope_index = pd.MultiIndex.from_arrays(
            arrays, names=["atomic_number", "mass_number"]
        )
    input_isotopes = pd.DataFrame(columns=np.arange(2),index = isotope_index)
    input_isotopes[0] = 0.01
    input_isotopes[1] = 0.02
    pd.testing.assert_frame_equal(csvy.raw_isotope_abundance,input_isotopes,check_names = True)
    

#test decay
def test_decay_dataframe(csvy):
    decay_index = pd.Index([1, 2, 26, 27, 28],name = 'atomic_number')
    ref_abund = pd.DataFrame([[0.00,0.4],[0.99,0.58],[0.00073564,0.00147129],[0.0076531,00.015306],[0.0016112,0.003222]], index = decay_index)

    pd.testing.assert_frame_equal(csvy.abundance,ref_abund,check_names = True, check_less_precise = True)
