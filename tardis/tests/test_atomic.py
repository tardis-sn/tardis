from tardis import atomic
from numpy import testing
import pytest
import os

def test_atomic_h5_readin():
    data = atomic.read_basic_atom_data(atomic.default_atom_h5_path)
    assert data['atomic_number'][13] == 14
    assert data['symbol'][13] == "Si"
    si_mass = data['mass'][13]
    testing.assert_almost_equal(si_mass, 28.085, decimal=4)
    pass

def test_ionization_h5_readin():
    data = atomic.read_ionization_data(atomic.default_atom_h5_path)
    hi_ionization = data['ionization_energy'][0]
    testing.assert_almost_equal(hi_ionization, 13.59844, decimal=4)
    
def test_lines_h5_readin():
    data=atomic.read_lines_data(atomic.default_atom_h5_path)
    assert data['line_id'][7] == 20
    assert data['wavelength'][7] == 71.533
    assert data['atomic_number'][7] == 14
    assert data['ion_number'][7] == 5
    assert data['level_number_lower'][7] == 1.0
    assert data['level_number_upper'][7] == 32.0
    si_frequency_ul = data['f_ul'][7]
    si_frequency_lu = data['f_lu'][7]
    testing.assert_almost_equal(si_frequency_ul,0.0772574,decimal=7)
    testing.assert_almost_equal(si_frequency_lu,0.1545148,decimal=7)

def test_collision_h5_readin():
    with pytest.raises(ValueError) as excinfo:
	atomic.read_collision_data(None)
    assert 'fname can not' in str(excinfo.value)
    test='../doesnt_exist.h5'
    with pytest.raises(IOError) as excinfo:
        atomic.read_collision_data(test)
    assert 'HDF5 File' in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
	atomic.read_collision_data(atomic.default_atom_h5_path)
    assert 'NLTE' in str(excinfo.value)
    path=os.path.join(os.path.dirname(__file__), 'data', 'chianti_he_db.h5')
    collision_data, temperature=atomic.read_collision_data(path)
    assert collision_data[0][0] == 2
    assert collision_data[0][1] == 0
    assert collision_data[0][2] == 18
    assert collision_data[0][3] == 2
    assert collision_data[0][4] == 1.0
    delta_e = collision_data[0][5]
    testing.assert_almost_equal(delta_e,35484.25114336,decimal=7)
    
def test_atom_levels():
    atom_data = atomic.AtomData.from_hdf5(atomic.default_atom_h5_path)
    with pytest.raises(Exception):
        raise Exception('test the atom_data thoroughly')

def test_atomic_symbol():
    assert atomic.atomic_number2symbol[14] == 'Si'

def test_atomic_symbol_reverse():
    assert atomic.symbol2atomic_number['Si'] == 14

@pytest.mark.skipif(not pytest.config.getvalue("atomic-dataset"),
                    reason='--atomic_database was not specified')
def test_atomic_reprepare():
    atom_data_filename = os.path.expanduser(os.path.expandvars(
        pytest.config.getvalue('atomic-dataset')))
    assert os.path.exists(atom_data_filename), ("{0} atomic datafiles "
                                                         "does not seem to "
                                                         "exist".format(
        atom_data_filename))
    atom_data = atomic.AtomData.from_hdf5(atom_data_filename)
    atom_data.prepare_atom_data([14])
    assert len(atom_data.lines) > 0
    atom_data.prepare_atom_data([20])
    assert len(atom_data.lines) > 0

