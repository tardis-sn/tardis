from tardis import atomic
from numpy import testing
import pytest
import h5py

@pytest.yield_fixture(scope="module")
def lines_dataset(atomic_data_fname):
    """ The fixture returns the dataset containing lines data from the provided dataset"""
    with h5py.File(atomic_data_fname, 'r') as h5_file:
        yield h5_file['lines_data']

@pytest.fixture
def atom_data_f_dataset_not_prepared(atomic_data_fname):
    """The fixture returns the AtomData instance that contains atomic data from the provided dataset.
       Does NOT call AtomData.prepare_atom_data method on the instance."""
    atom_data = atomic.AtomData.from_hdf5(atomic_data_fname)
    return atom_data

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


def test_atom_levels():
    atom_data = atomic.AtomData.from_hdf5(atomic.default_atom_h5_path)
    with pytest.raises(Exception):
        raise Exception('test the atom_data thoroughly')

def test_atomic_symbol():
    assert atomic.atomic_number2symbol[14] == 'Si'

def test_atomic_symbol_reverse():
    assert atomic.symbol2atomic_number['Si'] == 14


param_atom_num = [[14], [20], [14, 20]]

@pytest.mark.parametrize("selected_atomic_numbers", param_atom_num,
                         ids = ["atom_num: {}".format(_) for _ in param_atom_num])
def test_prepare_atom_data_set_lines(selected_atomic_numbers, atom_data_f_dataset_not_prepared, lines_dataset):
    """ Test that lines data is prepared in accordance with the selected atomic numbers
        Uses fixtures:
        --------
        atom_data_from_dataset : '~tardis.atomic.AtomData' instance containing data from the provided atomic dataset
        lines_dataset          :  HDF5 dataset "lines_data"

    """
    atom_data = atom_data_f_dataset_not_prepared
    atom_data.prepare_atom_data(selected_atomic_numbers)
    num_of_lines = 0
    # Go through the dataset and count the number of lines that should be selected
    for atom_num in lines_dataset['atomic_number']:
        if atom_num in selected_atomic_numbers:
            num_of_lines += 1

    assert len(atom_data.lines) == num_of_lines

param_atom_num_max_ion_num = [([14], 1), ([14, 20], 3)]

@pytest.mark.parametrize("selected_atomic_numbers, max_ion_number", param_atom_num_max_ion_num,
                         ids = ["atom_num: {}, max_ion_num: {}".format(*_) for _ in param_atom_num_max_ion_num])
def test_prepare_atom_data_set_lines_w_max_ion_number(selected_atomic_numbers, max_ion_number,
                                                      atom_data_f_dataset_not_prepared, lines_dataset):
    """ Test that lines data is prepared in accordance with the selected atomic numbers and the maximum ion number."""
    atom_data = atom_data_f_dataset_not_prepared
    atom_data.prepare_atom_data(selected_atomic_numbers, max_ion_number=max_ion_number)
    num_of_lines = 0
    for atom_num, ion_num in lines_dataset['atomic_number', 'ion_number']:
        if atom_num in selected_atomic_numbers and ion_num <= max_ion_number:
            num_of_lines += 1

    assert len(atom_data.lines) == num_of_lines

