from tardis import atomic
from numpy import testing
from astropy import units


def test_atomic_h5_readin():
    data = atomic.read_basic_atom_data()
    assert data['atomic_number'][13] == 14
    assert data['symbol'][13] == "Si"
    si_mass = units.Unit('g').to('u', data['mass'][13])
    testing.assert_almost_equal(si_mass, 28.085, decimal=4)
    pass

def test_ionization_h5_readin():
    data = atomic.read_ionization_data()
    HI_ionization = units.Unit('erg').to('eV', data['ionization_energy'][0])
    testing.assert_almost_equal(HI_ionization, 13.59844, decimal=4)
    pass

def test_atom_levels():
    atom_data = atomic.AtomData.from_hdf5()
    raise Exception('test the atom_data thoroughly')
    pass