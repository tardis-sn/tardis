from numpy import testing
import numpy as np
import pandas as pd
import pytest
from tardis import plasma_array, atomic
from tardis.util import intensity_black_body
from tardis.io.config_reader import ConfigurationNameSpace
import os
import tardis
from astropy import units as u, constants as const

data_path = os.path.join(tardis.__path__[0], 'tests', 'data')
helium_test_db = os.path.join(data_path, 'chianti_he_db.h5')

ion_populations=pd.DataFrame([[2,0,1.0],[2,1,1.0],[2,2,1.0]],
                             columns=['atomic_number', 'ion_number', 0])
ion_populations = ion_populations.set_index(['atomic_number', 'ion_number'])
from astropy import constants as const, units as u
from numpy.testing import assert_allclose


def calculate_lte_level_populations(atom_data, ion_number, temperature):
    #self.atom_data.levels["g"].ix[(2,ion_number)].values * np.exp(- self.atom_data.levels["energy"].ix[(2,ion_number)].values * u.erg / const.k_B / self.plasma.t_rads / u.K).value
    g = atom_data.levels["g"].ix[(2,ion_number)].values
    energy = atom_data.levels["energy"].ix[(2,ion_number)].values  * u.erg
    temperature = temperature * u.K
    level_populations = (g * np.exp(-energy / (temperature * const.k_B))).value
    return level_populations

def pytest_generate_tests(metafunc):
    # called once per each test function
    funcarglist = metafunc.cls.params[metafunc.function.__name__]
    argnames = list(funcarglist[0])
    metafunc.parametrize(argnames, [[funcargs[name] for name in argnames]
            for funcargs in funcarglist])

class TestNLTELTEApproximation(object):

    params = {"test_He_ltelevelpops" : [dict(ion_number = 0),
                                        dict(ion_number = 1)]}

    def setup(self):
        self.nlte_species=[(2,0),(2,1)]
        self.nlte_config = ConfigurationNameSpace({'species':
                                                             self.nlte_species})
        self.atom_data = atomic.AtomData.from_hdf5(helium_test_db)

        self.plasma = plasma_array.BasePlasmaArray.from_abundance(
            {'He':1.0}, 1e-15 * u.Unit('g/cm3'), self.atom_data, 10 * u.day,
            nlte_config=self.nlte_config)
        self.plasma.j_blues = pd.DataFrame(intensity_black_body(
            self.atom_data.lines.nu.values[np.newaxis].T, np.array([10000.])))
        self.plasma.tau_sobolevs = pd.DataFrame(np.zeros_like(
            self.plasma.j_blues))
        self.plasma.t_rads=np.array([10000.])
        self.plasma.t_electrons=np.array([10000.])
        self.plasma.ws=np.array([1.0])
        self.plasma.electron_densities=pd.Series([1.e9])
        self.plasma.ion_populations = ion_populations
        self.plasma.calculate_nlte_level_populations()

    def test_He_ltelevelpops(self, ion_number):
        lte_pops = calculate_lte_level_populations(self.atom_data, ion_number,
                                                   self.plasma.t_rads)


        g_ground_level = self.atom_data.levels["g"].ix[(2,ion_number)][0]
        np.testing.assert_allclose(lte_pops,
                                   g_ground_level *
                                   self.plasma.level_populations[0].
                                   ix[(2,ion_number)].values)

class TestNLTE(object):

    params = {"test_He_dilutelevelpops" : [dict(dummy = 0) ],
              "test_He_dilutelevelpops_isnotLTE" : [dict(ion_number = 0),
                                                    dict(ion_number = 1)]}

    def setup(self):
        self.nlte_species=[(2,0),(2,1)]
        self.nlte_config = ConfigurationNameSpace({'species':self.nlte_species})
        self.atom_data = atomic.AtomData.from_hdf5(helium_test_db)
        self.plasma = plasma_array.BasePlasmaArray.from_abundance(
            {'He':1.0}, 1e-15*u.Unit('g/cm3'), self.atom_data, 10 * u.day,
            nlte_config=self.nlte_config)
        self.plasma.j_blues = 0.5 * pd.DataFrame(intensity_black_body(self.atom_data.lines.nu.values[np.newaxis].T, np.array([10000.])))
        self.plasma.tau_sobolevs = pd.DataFrame(np.zeros_like(self.plasma.j_blues))
        self.plasma.t_rads=np.array([10000.])
        self.plasma.t_electrons=np.array([9000.])
        self.plasma.ws=np.array([0.5])
        self.plasma.electron_densities=pd.Series([1.e9])
        self.plasma.ion_populations = ion_populations
        self.plasma.calculate_nlte_level_populations()

    def test_He_dilutelevelpops(self, dummy):
        ref_pops = pd.read_hdf(os.path.join(data_path,
                                            'He_nlte_pops.h5'), 'He_level_pops')
        np.testing.assert_allclose(self.plasma.level_populations.values, ref_pops.values)

    def test_He_dilutelevelpops_isnotLTE(self, ion_number):
        lte_pops = self.atom_data.levels["g"].ix[(2,ion_number)].values * np.exp(- self.atom_data.levels["energy"].ix[(2,ion_number)].values * u.erg / const.k_B / self.plasma.t_rads / u.K).value
        assert not np.allclose(lte_pops, self.atom_data.levels["g"].ix[(2,ion_number)][0]*self.plasma.level_populations[0].ix[(2,ion_number)].values, atol=0)


