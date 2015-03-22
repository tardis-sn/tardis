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

def test_lines_h5_readin()
    data = atomic.read_lines_data(atomic.default_atom_h5_path)
    assert data['line_id'][8] == 21
    assert data['wavelength'][8] == 71.5610
    B_ul = data['B_ul'][8]
    A_ul = data['A_ul'][8]
    testing.assert_almost_equal(A_ul,34179006146.602795,decimal=4)
    testing.assert_almost_equal(B_ul,31526948.884054925,decimal=4)
    pass


def test_levels_h5_readin():
    data = atomic.read_levels_data(atomic.default_atom_h5_path)
    assert data['atomic_number'][4] == 14
    assert data['ion_number'][4] == 0
    assert data['level_number'][4] == 4
    si_energy = data['energy'][4]
    testing.assert_almost_equal(si_energy, 1.90865, decimal=4)
    assert data['g'][4] == 1
    assert data['metastable'][4] == False


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

