from astropy import constants as const, units as u
import os
import pandas as pd
import tardis
from tardis import atomic
import warnings
from tardis.plasma.standard_plasmas import LegacyPlasmaArray
from tardis.io.util import parse_abundance_dict_to_dataframe
# from numpy.testing import assert_allclose
data_path = os.path.join(tardis.__path__[0], 'tests', 'data')
helium_test_db = os.path.join(data_path, 'chianti_he_db.h5')


class TestNebularPlasma(object):

    def setup(self):
        atom_data = atomic.AtomData.from_hdf5(helium_test_db)
        density = 1e-15 * u.Unit('g/cm3')
        abundance = parse_abundance_dict_to_dataframe({'He': 1.0})
        abundance = pd.DataFrame({0: abundance})
        atom_data.prepare_atom_data([2])  # FIXME: hardcoded, bad
        number_densities = abundance * density.to('g/cm^3').value
        number_densities = number_densities.div(
            atom_data.atom_data.mass.ix[number_densities.index], axis=0)
        self.plasma = LegacyPlasmaArray(
                number_densities,
                atomic_data=atom_data, time_explosion=10 * u.day,
                ionization_mode='nebular')

#        self.plasma = plasma_array.BasePlasmaArray.from_abundance(
#            {'He':1.0}, 1e-15*u.Unit('g/cm3'), atom_data, 10 * u.day,
#            ionization_mode='nebular', excitation_mode='dilute-lte')

    def test_high_temperature(self):
        with warnings.catch_warnings(record=True) as w:
            self.plasma.update_radiationfield(
                    t_rad=[100000.], ws=[0.5], j_blues=None, nlte_config=None)
        assert str(w[0].message).startswith(
                't_rads outside of zeta factor interpolation')
