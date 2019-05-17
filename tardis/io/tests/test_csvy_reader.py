from tardis.io.parsers import csvy
import pytest
import os
from astropy import units as u
import numpy.testing as npt

test_data_directory = os.path.dirname(__file__)

def data_path(filename):
    data_dir = os.path.dirname(__file__)
    return os.path.join(data_dir, 'data', filename)

def test_csvy_finds_csv_first_line():
    yaml_dict, csv = csvy.load_csvy(data_path('csvy_full.csvy'))
    assert csv['velocity'][0] == 10000

def test_csv_colnames_equiv_datatype_fields():
    yaml_dict, csv = csvy.load_csvy(data_path('csvy_full.csvy'))
    datatype_names = [od['name'] for od in yaml_dict['datatype']['fields']]
    for key in csv.columns:
        assert key in datatype_names
    for name in datatype_names:
        assert name in csv.columns

def test_csvy_nocsv_data_is_none():
    yaml_dict, csv = csvy.load_csvy(data_path('csvy_nocsv.csvy'))
    assert csv is None
