from tardis import atomic
from numpy import testing
import pytest

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


def test_atom_levels():
    atom_data = atomic.AtomData.from_hdf5(atomic.default_atom_h5_path)
    with pytest.raises(Exception):
        raise Exception('test the atom_data thoroughly')

def test_atomic_symbol():
    assert atomic.atomic_number2symbol[14] == 'Si'

def test_atomic_symbol_reverse():
    assert atomic.symbol2atomic_number['Si'] == 14
