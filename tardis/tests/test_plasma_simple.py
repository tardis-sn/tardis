from astropy import constants as const, units as u
import os
import tardis
from tardis import plasma_array, atomic
import pytest
#from numpy.testing import assert_allclose
data_path = os.path.join(tardis.__path__[0], 'tests', 'data')
helium_test_db = os.path.join(data_path, 'chianti_he_db.h5')


class TestNebularPlasma(object):

    def setup(self):
        atom_data = atomic.AtomData.from_hdf5(helium_test_db)
        self.plasma = plasma_array.BasePlasmaArray.from_abundance(
            {'He':1.0}, 1e-15*u.Unit('g/cm3'), atom_data, 10 * u.day,
            ionization_mode='nebular', excitation_mode='dilute-lte')

    def test_high_temperature(self):
        with pytest.raises(ValueError) as excinfo:
            self.plasma.update_radiationfield([100000.], [1.])

        assert str(excinfo.value).startswith('t_rads outside of zeta '
                                                'factor interpolation')