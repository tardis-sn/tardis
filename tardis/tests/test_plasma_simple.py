from astropy import constants as const, units as u
import os
import pandas as pd
import tardis
from tardis.io.atomic import AtomData
import warnings
import pytest
from tardis.io.util import parse_abundance_dict_to_dataframe


# FIXME
@pytest.mark.skipif(
        True,
        reason="This test needs a rewrite since LegacyPlasmaArray does no"
        "longer exist.")
class TestNebularPlasma(object):

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        atom_data_filename = os.path.expanduser(os.path.expandvars(
            pytest.config.getvalue('atomic-dataset')))
        assert os.path.exists(atom_data_filename), \
            "{0} atomic datafiles does not seem to exist".format(atom_data_filename)

        atom_data = AtomData.from_hdf(atom_data_filename)
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
