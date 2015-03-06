from tardis import atomic
from numpy import testing, float64
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

def test_levels_h5_readin():
    data = atomic.read_levels_data(atomic.default_atom_h5_path)
    assert data['atomic_number'][4] == 14
    assert data['ion_number'][4] == 0
    assert data['level_number'][4] == 4
    si_energy = data['energy'][4]
    testing.assert_almost_equal(si_energy, 1.90865, decimal=4)
    assert data['g'][4] == 1
    assert data['metastable'][4] == False

def test_lines_h5_readin():
    data = atomic.read_lines_data(atomic.default_atom_h5_path)
    assert data['line_id'][5] == 18
    si_wavelength = data['wavelength'][5]
    testing.assert_almost_equal(si_wavelength, 71.273)
    assert data['atomic_number'][5] == 14
    assert data['ion_number'][5] == 5
    si_f_ul = data['f_ul'][5]
    testing.assert_almost_equal(si_f_ul, 0.00075, decimal=4)
    si_f_lu = data['f_lu'][5]
    testing.assert_almost_equal(si_f_lu, 0.00075, decimal=4)
    si_level_number_lower = data['level_number_lower'][5]
    testing.assert_almost_equal(si_level_number_lower, 0.0)
    si_level_number_upper = data['level_number_upper'][5]
    testing.assert_almost_equal(si_level_number_upper, 32.0)
    si_nu = data['nu'][5]
    testing.assert_almost_equal(si_nu, float64(4.2062556367769e+16))
    si_B_lu = data['B_lu'][5]
    testing.assert_almost_equal(si_B_lu, 903442.64270, decimal=4)
    si_B_ul = data['B_ul'][5]
    testing.assert_almost_equal(si_B_ul, 903442.64270, decimal=4)
    si_A_ul = data['A_ul'][5]
    testing.assert_almost_equal(si_A_ul, 991361747.81935, decimal=4)


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

