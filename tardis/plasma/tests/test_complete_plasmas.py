import os

import pandas as pd
import pytest
from pandas.util import testing as pdt

from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation


class BasePlasmaTest:
    #Class defining all common tests for different setups of Plasma
    #This can then be inherited for different Plasma setup
    @classmethod
    def setup(cls):
        pass

    def read_hdf_attr(self, attr):
        return pd.read_hdf(self.reference_file_path, os.path.join('plasma', attr))

    #GENERAL PROPERTIES
    def test_beta_rad(self):
        pdt.assert_almost_equal(self.plasma.beta_rad,
                                self.read_hdf_attr('beta_rad').values)

    def test_g_electron(self):
        pdt.assert_almost_equal(self.plasma.g_electron,
                                self.read_hdf_attr('g_electron').values)

    def test_number_density(self):
        pdt.assert_almost_equal(
            self.plasma.number_density, self.read_hdf_attr('number_density'))

    def test_electron_temperature(self):
        pdt.assert_almost_equal(self.plasma.t_electrons,
                                self.read_hdf_attr('t_electrons').values)

    #PARTITION PROPERTIES
    def test_level_boltzmann_factor(self):
        pdt.assert_almost_equal(
            self.plasma.level_boltzmann_factor, self.read_hdf_attr('level_boltzmann_factor'))

    def test_lte_partition_function(self):
        pdt.assert_almost_equal(
            self.plasma.partition_function, self.read_hdf_attr('partition_function'))

    #ATOMIC PROPERTIES
    def test_levels_property(self):
        pdt.assert_almost_equal(
            self.plasma.excitation_energy, self.read_hdf_attr('excitation_energy'))

    def test_lines_property(self):
        pdt.assert_almost_equal(self.plasma.lines, self.read_hdf_attr('lines'))

    def test_lines_lower_level_index_property(self):
        pdt.assert_almost_equal(self.plasma.lines_lower_level_index, self.read_hdf_attr(
            'lines_lower_level_index').values)

    def test_lines_upper_level_index_property(self):
        pdt.assert_almost_equal(self.plasma.lines_upper_level_index, self.read_hdf_attr(
            'lines_upper_level_index').values)

    def test_atomic_mass_property(self):
        pdt.assert_almost_equal(self.plasma.atomic_mass,
                                self.read_hdf_attr('atomic_mass'))

    def test_ionization_data_property(self):
        pdt.assert_almost_equal(
            self.plasma.ionization_data, self.read_hdf_attr('ionization_data'))

    #ION POPULATION PROPERTIES
    def test_phi(self):
        pdt.assert_almost_equal(self.plasma.phi, self.read_hdf_attr('phi'))

    def test_ion_number_density(self):
        pdt.assert_almost_equal(
            self.plasma.ion_number_density, self.read_hdf_attr('ion_number_density'))

    def test_electron_densities(self):
        pdt.assert_almost_equal(
            self.plasma.electron_densities, self.read_hdf_attr('electron_densities'))

    #LEVEL POPULATION PROPERTIES
    def test_level_number_density(self):
        pdt.assert_almost_equal(
            self.plasma.level_number_density, self.read_hdf_attr('level_number_density'))

    #RADIATIVE PROPERTIES
    def test_stimulated_emission_factor(self):
        pdt.assert_almost_equal(self.plasma.stimulated_emission_factor, self.read_hdf_attr(
            'stimulated_emission_factor').values)

    def test_tau_sobolev(self):
        pdt.assert_almost_equal(self.plasma.tau_sobolevs,
                                self.read_hdf_attr('tau_sobolevs'))

    def test_beta_sobolev(self):
        pdt.assert_almost_equal(self.plasma.beta_sobolev,
                                self.read_hdf_attr('beta_sobolev').values)


@pytest.mark.skipif(not pytest.config.getvalue("tardis-refdata"),
                    reason="Path to Tardis Reference files is not defined")
class TestLTEPlasma(BasePlasmaTest):

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(cls):
        tardis_ref_path = os.path.expanduser(
            os.path.expandvars(pytest.config.getvalue('tardis-refdata')))
        cls.reference_file_path = os.path.join(
            tardis_ref_path, 'plasma_reference', 'plasma_lte_reference.h5')

        cls.config_path = os.path.join(
            'tardis', 'plasma', 'tests', 'data', 'plasma_test_config_lte.yml')
        cls.config = Configuration.from_yaml(cls.config_path)
        cls.config['atom_data'] = os.path.join(
            tardis_ref_path, 'atom_data', 'kurucz_cd23_chianti_H_He.h5')
        cls.sim = Simulation.from_config(cls.config)
        cls.sim.run()
        cls.plasma = cls.sim.plasma


@pytest.mark.skipif(not pytest.config.getvalue("tardis-refdata"),
                    reason="Path to Tardis Reference files is not defined")
class TestNLTEPlasma(BasePlasmaTest):

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(cls):
        tardis_ref_path = os.path.expanduser(
            os.path.expandvars(pytest.config.getvalue('tardis-refdata')))
        cls.reference_file_path = os.path.join(
            tardis_ref_path, 'plasma_reference', 'plasma_nlte_reference.h5')

        cls.config_path = os.path.join(
            'tardis', 'plasma', 'tests', 'data', 'plasma_test_config_nlte.yml')
        cls.config = Configuration.from_yaml(cls.config_path)
        cls.config['atom_data'] = os.path.join(
            tardis_ref_path, 'atom_data', 'kurucz_cd23_chianti_H_He.h5')
        cls.sim = Simulation.from_config(cls.config)
        cls.sim.run()
        cls.plasma = cls.sim.plasma

        #Additional Tests for NLTE Plasma, apart from tests defined in BasePlasmaTest
        #GENERAL PROPERTIES
        def test_beta_electron(self):
            pdt.assert_almost_equal(self.plasma.beta_electron,
                                    self.read_hdf_attr('beta_electron').values)

        def test_selected_atoms(self):
            pdt.assert_almost_equal(
                self.plasma.selected_atoms.values, self.read_hdf_attr('selected_atoms').values)

        #ATOMIC PROPERTIES
        def test_zeta_data_property(self):
            pdt.assert_almost_equal(
                self.plasma.zeta_data.values, self.read_hdf_attr('zeta_data').values)

        #ION POPULATION PROPERTIES
        def test_radiation_field_correction(self):
            pdt.assert_almost_equal(
                self.plasma.delta, self.read_hdf_attr('delta'))
