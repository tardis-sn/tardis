__author__ = 'maryam'

#this is the test environment for atomic

from tardis import atomic
from numpy import testing
from astropy import units


def test_atomic_h5_readin():
    data = atomic.read_basic_atom_data()
    assert data['atomic_number'][13] == 14
    assert data['symbol'][13] == "Si"
    si_mass = units.Unit('g').to('u', data['mass'][13])
    testing.assert_almost_equal(si_mass, 28.085, decimal=4)


def test_ionization_h5_readin():
    data = atomic.read_ionization_data()
    HI_ionization = units.Unit('erg').to('eV', data['ionization_energy'][0])
    testing.assert_almost_equal(HI_ionization, 13.59844, decimal=4)


def test_atom_levels():
    atom_data = atomic.AtomData.from_hdf5()._levels
    assert atom_data[0]['level_number'] == atom_data[0]['energy']
    for i in range(len(atom_data)):
        assert atom_data[i][0] == 14
        assert atom_data[i][1] == 1
        assert atom_data[i][2] == i

def test_atom_lines():
    atom_data = atomic.AtomData.from_hdf5()
    wavelength = atom_data._lines['wavelength'][0]
    testing.assert_almost_equal(wavelength, 711.338999, decimal=4)
    assert atom_data._lines['atomic_number'][0] == 14
    assert atom_data._lines['ion_number'][0] == 1
    f_ul = atom_data._lines['f_ul'][0]
    testing.assert_almost_equal(f_ul, 0.018532, decimal=4)
    f_lu = atom_data._lines['f_lu'][0]
    testing.assert_almost_equal(f_lu, 0.037065, decimal=4)
    assert atom_data._lines['level_id_lower'][0] == 0
    assert atom_data._lines['level_id_upper'][0] == 113


