import os
import pytest
from numpy import testing
from tardis.io import atomic


@pytest.fixture(scope="module")
def default_atom_h5_path():
    return atomic.default_atom_h5_path


@pytest.fixture(scope="module")
def chianti_he_db_h5_path():
    return os.path.join(
            os.path.dirname(os.path.realpath(__file__)), 'data', 'chianti_he_db.h5')


def test_data_path():
    data_path = atomic.data_path('test')
    assert data_path.split('/')[-3:] == ['tardis', 'data', 'test']


def test_read_basic_atom_data(default_atom_h5_path):
    data = atomic.read_basic_atom_data(default_atom_h5_path)
    assert data['atomic_number'][13] == 14
    assert data['symbol'][13] == "Si"
    testing.assert_almost_equal(data['mass'][13], 28.085, decimal=4)


def test_read_ionization_data(default_atom_h5_path):
    data = atomic.read_ionization_data(default_atom_h5_path)
    assert data['atomic_number'][0] == 1
    assert data['ion_number'][0] == 1
    testing.assert_almost_equal(data['ionization_energy'][0], 13.59844, decimal=4)


def test_read_levels_data(default_atom_h5_path):
    data = atomic.read_levels_data(default_atom_h5_path)
    assert data['atomic_number'][4] == 14
    assert data['ion_number'][4] == 0
    assert data['level_number'][4] == 4
    testing.assert_almost_equal(data['energy'][4], 1.90865, decimal=4)
    assert data['g'][4] == 1
    assert data['metastable'][4] == False


def test_read_lines_data(default_atom_h5_path):
    data = atomic.read_lines_data(default_atom_h5_path)
    assert data['line_id'][0] == 8
    assert data['atomic_number'][0] == 14
    assert data['ion_number'][0] == 5
    testing.assert_almost_equal(data['wavelength'][0], 66.772, decimal=4)
    testing.assert_almost_equal(data['f_ul'][0], 0.02703, decimal=4)
    testing.assert_almost_equal(data['f_lu'][0], 0.04054, decimal=4)
    assert data['level_number_lower'][0] == 0.0
    assert data['level_number_upper'][0] == 36.0


def test_read_synpp_refs(chianti_he_db_h5_path):
    data = atomic.read_synpp_refs(chianti_he_db_h5_path)
    assert data['atomic_number'][0] == 1
    assert data['ion_number'][0] == 0
    testing.assert_almost_equal(data['wavelength'][0], 6562.7973633, decimal=4)
    assert data['line_id'][0] == 564995


def test_read_zeta_data(default_atom_h5_path, chianti_he_db_h5_path):
    data = atomic.read_zeta_data(chianti_he_db_h5_path)
    testing.assert_almost_equal(data[2000][1][1], 0.339, decimal=4)
    testing.assert_almost_equal(data[2000][1][2], 0.000, decimal=4)

    with pytest.raises(ValueError):
        atomic.read_zeta_data(None)

    with pytest.raises(IOError):
        atomic.read_zeta_data('fakepath')

    with pytest.raises(ValueError):
        atomic.read_zeta_data(default_atom_h5_path)


def test_read_collision_data(default_atom_h5_path, chianti_he_db_h5_path):
    data = atomic.read_collision_data(chianti_he_db_h5_path)
    assert data[0]['atomic_number'][0] == 2
    assert data[0]['ion_number'][0] == 0
    assert data[0]['level_number_upper'][0] == 18
    assert data[0]['level_number_lower'][0] == 2
    assert data[0]['g_ratio'][0] == 1.0
    testing.assert_almost_equal(data[0]['delta_e'][0], 35484.251143, decimal=4)
    assert data[1][0] == 2000.0
    assert data[1][1] == 4000.0

    with pytest.raises(ValueError):
        atomic.read_zeta_data(None)

    with pytest.raises(IOError):
        atomic.read_zeta_data('fakepath')

    with pytest.raises(ValueError):
        atomic.read_zeta_data(default_atom_h5_path)


def test_read_macro_atom_data(default_atom_h5_path, chianti_he_db_h5_path):
    data = atomic.read_macro_atom_data(chianti_he_db_h5_path)
    assert data[0]['atomic_number'][0] == 2
    assert data[0]['ion_number'][0] == 0
    assert data[0]['source_level_number'][0] == 0.0
    assert data[0]['destination_level_number'][0] == 48.0
    assert data[0]['transition_type'][0] == 1
    assert data[0]['transition_probability'][0] == 0.0
    assert data[0]['transition_line_id'][0] == 564957

    assert data[1]['count_down'][0] == 0
    assert data[1]['count_up'][0] == 7
    assert data[1]['count_total'][0] == 7

    with pytest.raises(ValueError):
        atomic.read_macro_atom_data(None)

    with pytest.raises(IOError):
        atomic.read_macro_atom_data('fakepath')

    with pytest.raises(ValueError):
        atomic.read_macro_atom_data(default_atom_h5_path)


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
    # Fix for new behavior of prepare_atom_data
    # Consider doing only one prepare_atom_data and check
    # len(atom_data.lines) == N where N is known
    atom_data = atomic.AtomData.from_hdf5(atom_data_filename)
    atom_data.prepare_atom_data([20])
    assert len(atom_data.lines) > 0

