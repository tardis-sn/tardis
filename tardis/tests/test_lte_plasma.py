from tardis import plasma_array as plasma, atomic
from astropy import constants
import numpy as np
import pytest

atom_data = atomic.AtomData.from_hdf5(atomic.default_atom_h5_path)
atom_data.prepare_atom_data(selected_atomic_numbers=[14])

pytestmark = pytest.mark.skipif(True, reason='to be implemented')


class TestNormalLTEPlasma:

    compare_part_func = np.array([11.406201367482032, 5.866632552894803, 1.0044215520812598, 2.0002017142942163,
                                        1.0, 4.961567712516646])

    compare_phis = np.array([9.679214946588846e+16, 2390247610210.469, 63428.65716485618, 0.021444877147797966,
                             1.0749561507478891e-62])

    compare_ion_populations = np.array([0.15945011675589407, 3717655.969887028, 2140505324.5794258, 32704.389884367,
                                        1.6894052370987268e-07, 0.0])

    compare_level_populations_14_1 = np.array([1267390.086686989, 2432159.357552647, 2673.1427419153792,
                                               5263.602861186024, 7698.769675188349])

    def setup(self):
        self.plasma = plasma.LTEPlasma.from_abundance(10000, {'Si':1.0}, 1e-13, atom_data, 10*86400)

    def test_beta_rad(self):
        assert self.plasma.beta_rad == 1 / (10000 * constants.k_B.cgs.value)

    def test_t_electron(self):
        assert self.plasma.t_electron == 0.9 * self.plasma.t_rad

    def test_saha_calculation_method(self):
        assert self.plasma.calculate_saha == self.plasma.calculate_saha_lte

    def test_partition_function_calculation(self):
        assert np.all(self.plasma.partition_functions.values == self.compare_part_func)

    def test_phis_calculation(self):
        self.calculated_phis = self.plasma.calculate_saha()

        assert np.all(self.calculated_phis.values == self.compare_phis)

    def test_ionization_balance_calculation(self):
        assert np.all(self.plasma.ion_populations.values == self.compare_ion_populations)

    def test_electron_density(self):
        assert self.plasma.electron_density == 4284826418.2983923

    def test_level_populations(self):
        assert np.all(self.plasma.level_populations.ix[14, 1].values[:5] == self.compare_level_populations_14_1)

    def test_tau_sobolev(self):
        #silicon line 14, 1 , wl = 6347.105178 level_lower = 7 level_upper = 12
        wavelength_id = 565376
        assert self.plasma.tau_sobolevs[self.plasma.atom_data.lines.index == wavelength_id][0] == 101.06456251838634


    def test_population_inversion(self):
        self.plasma.level_populations.ix[14, 1, 12] = 1.1 * self.plasma.level_populations.ix[14, 1, 7]
        with pytest.raises(plasma.PopulationInversionException):
            self.plasma.calculate_tau_sobolev()


