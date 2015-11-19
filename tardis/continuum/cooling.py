from tardis.continuum.base import ContinuumProcess
import numpy as np


class CoolingRates(ContinuumProcess):
    def __init__(self, input_data, coll_excitation_cooling_rate, coll_ionization_cooling_rate, fb_cooling_rate):
        super(CoolingRates, self).__init__(input_data)

        self.ff_cooling = \
            self._calculate_ff_cooling_rate()
        self.collisional_excitation_cooling = coll_excitation_cooling_rate
        self.collisional_ionization_cooling = coll_ionization_cooling_rate
        self.fb_cooling = fb_cooling_rate
        self._set_total_cooling_rates()
        self._set_cooling_probabilities()
        self._set_individual_cooling_probabilities()
        self._prepare_montecarlo_data()

    def _calculate_ff_cooling_rate(self):
        C0 = self.input.C0_ff
        # TODO: value for Gaunt factor (Lucy: = 1; Osterbrock recommendation for nebular conditions: = 1.3 )
        factor = self.ion_number_density.mul(np.square(self.input.ion_charges), axis=0).sum().values
        cooling_rate = C0 * self.electron_densities * np.sqrt(self.t_electrons) * factor
        return cooling_rate

    def _set_total_cooling_rates(self):
        self.collisional_excitation_cooling_total = self.collisional_excitation_cooling.sum().values
        self.collisional_ionization_cooling_total = self.collisional_ionization_cooling.sum().values
        self.fb_cooling_total = self.fb_cooling.sum().values

    def _set_cooling_probabilities(self):
        total_cooling_rate = self.fb_cooling_total + self.ff_cooling + self.collisional_excitation_cooling_total + \
                             self.collisional_ionization_cooling_total
        self.fb_cooling_prob = self.fb_cooling_total / total_cooling_rate
        self.ff_cooling_prob = self.ff_cooling / total_cooling_rate
        self.coll_ion_cooling_prob = self.collisional_ionization_cooling_total / total_cooling_rate
        self.coll_exc_cooling_prob = self.collisional_excitation_cooling_total / total_cooling_rate

    def _set_individual_cooling_probabilities(self):
        self.fb_cooling_prob_individual = self.fb_cooling.divide(self.fb_cooling_total, axis=1)
        self.coll_exc_cooling_prob_individual = self.collisional_excitation_cooling.divide(
            self.collisional_excitation_cooling_total, axis=1)
        self.coll_ion_cooling_prob_individual = self.collisional_ionization_cooling.divide(
            self.collisional_ionization_cooling_total, axis=1)

    def _prepare_montecarlo_data(self):
        self.fb_cooling_prob_array = self._get_contiguous_array(self.fb_cooling_prob_individual)
        self.coll_exc_cooling_prob_array = self._get_contiguous_array(self.coll_exc_cooling_prob_individual)
        self.coll_ion_cooling_prob_array = self._get_contiguous_array(self.coll_ion_cooling_prob_individual)

    def _get_contiguous_array(self, dataframe):
        return np.ascontiguousarray(dataframe.values.transpose())

    @property
    def fb_cooling_prob_nd(self):
        return self.fb_cooling_prob_array.shape[1]

    @property
    def coll_exc_cooling_prob_nd(self):
        return self.coll_exc_cooling_prob_array.shape[1]

    @property
    def coll_ion_cooling_prob_nd(self):
        return self.coll_ion_cooling_prob_array.shape[1]