from collections import namedtuple

import numpy as np

from tardis.iip_plasma.continuum.base import ContinuumProcess

CoolingData = namedtuple(
    "CoolingData",
    [
        "cooling_probability",
        "probabilities_array",
        "prob_array_nd",
        "references",
    ],
)


class CoolingRates(ContinuumProcess):
    def __init__(self, input_data, **kwargs):
        super(CoolingRates, self).__init__(input_data)
        self.total_cooling_rate = np.zeros(self.no_of_shells)
        self.cooling_processes = kwargs.keys()
        self._set_cooling_rates(**kwargs)

    def _set_cooling_rates(self, **kwargs):
        total_cooling_rates_dict = {}
        for name, cooling_process in kwargs.iteritems():
            cooling_rate_of_process = cooling_process.cooling_rate
            total_cooling_rate_of_process = self._get_total_cooling_rate(
                cooling_process
            )
            total_cooling_rates_dict[name] = total_cooling_rate_of_process
            individual_cooling_probabilities = (
                self._get_individual_cooling_probabilities(
                    cooling_rate_of_process,
                    total_cooling_rate_of_process,
                    name=cooling_process.name,
                )
            )

            setattr(self, name + "cooling_rate", cooling_rate_of_process)
            setattr(self, name + "_total", total_cooling_rate_of_process)
            if individual_cooling_probabilities is not None:
                individual_cooling_probabilities_array = (
                    self._get_contiguous_array(individual_cooling_probabilities)
                )
                prob_array_nd = individual_cooling_probabilities_array.shape[1]
                references = individual_cooling_probabilities.index.values
                setattr(
                    self,
                    name + "_individual_probabilities",
                    individual_cooling_probabilities,
                )
                setattr(
                    self,
                    name,
                    CoolingData(
                        np.zeros(self.no_of_shells),
                        individual_cooling_probabilities_array,
                        prob_array_nd,
                        references,
                    ),
                )

            self.total_cooling_rate += total_cooling_rate_of_process

        self._set_cooling_probabilities(total_cooling_rates_dict)

    def _get_total_cooling_rate(self, cooling_process):
        if cooling_process.name == "free_free":
            return cooling_process.cooling_rate
        return cooling_process.cooling_rate.sum().values

    def _set_cooling_probabilities(self, total_cooling_rates_dict):
        for (
            name,
            total_cooling_rate_of_process,
        ) in total_cooling_rates_dict.iteritems():
            cooling_probability_of_process = (
                total_cooling_rate_of_process / self.total_cooling_rate
            )
            setattr(self, name + "_probability", cooling_probability_of_process)
            if not name == "free_free":
                process = getattr(self, name)
                process.cooling_probability[:] = cooling_probability_of_process

    def _get_individual_cooling_probabilities(
        self, cooling_rate_of_process, total_cooling_rate_of_process, name
    ):
        if not name == "free_free":
            return cooling_rate_of_process.divide(
                total_cooling_rate_of_process, axis=1
            )
        return None


class CoolingRatesLegacy(ContinuumProcess):
    def __init__(
        self,
        input_data,
        coll_excitation_cooling_rate,
        coll_ionization_cooling_rate,
        fb_cooling_rate,
    ):
        super(CoolingRatesLegacy, self).__init__(input_data)

        self.ff_cooling = self._calculate_ff_cooling_rate()
        self.collisional_excitation_cooling = coll_excitation_cooling_rate
        self.collisional_ionization_cooling = coll_ionization_cooling_rate
        self.fb_cooling = fb_cooling_rate
        self._set_total_cooling_rates()
        self._set_cooling_probabilities()
        self._set_individual_cooling_probabilities()
        self._set_montecarlo_data()

    def _calculate_ff_cooling_rate(self):
        C0 = self.input.C0_ff
        # TODO: value for Gaunt factor (Lucy: = 1; Osterbrock recommendation for nebular conditions: = 1.3 )
        factor = (
            self.ion_number_density.mul(
                np.square(self.input.ion_charges), axis=0
            )
            .sum()
            .values
        )
        cooling_rate = (
            C0 * self.electron_densities * np.sqrt(self.t_electrons) * factor
        )
        return cooling_rate

    def _set_total_cooling_rates(self):
        self.collisional_excitation_cooling_total = (
            self.collisional_excitation_cooling.sum().values
        )
        self.collisional_ionization_cooling_total = (
            self.collisional_ionization_cooling.sum().values
        )
        self.fb_cooling_total = self.fb_cooling.sum().values

    def _set_cooling_probabilities(self):
        total_cooling_rate = (
            self.fb_cooling_total
            + self.ff_cooling
            + self.collisional_excitation_cooling_total
            + self.collisional_ionization_cooling_total
        )
        self.fb_cooling_prob = self.fb_cooling_total / total_cooling_rate
        self.ff_cooling_prob = self.ff_cooling / total_cooling_rate
        self.coll_ion_cooling_prob = (
            self.collisional_ionization_cooling_total / total_cooling_rate
        )
        self.coll_exc_cooling_prob = (
            self.collisional_excitation_cooling_total / total_cooling_rate
        )

    def _set_individual_cooling_probabilities(self):
        self.fb_cooling_prob_individual = self.fb_cooling.divide(
            self.fb_cooling_total, axis=1
        )
        self.coll_exc_cooling_prob_individual = (
            self.collisional_excitation_cooling.divide(
                self.collisional_excitation_cooling_total, axis=1
            )
        )
        self.coll_ion_cooling_prob_individual = (
            self.collisional_ionization_cooling.divide(
                self.collisional_ionization_cooling_total, axis=1
            )
        )

    def _set_montecarlo_data(self):
        self.fb_cooling_prob_array = self._get_contiguous_array(
            self.fb_cooling_prob_individual
        )
        self.coll_exc_cooling_prob_array = self._get_contiguous_array(
            self.coll_exc_cooling_prob_individual
        )
        self.coll_ion_cooling_prob_array = self._get_contiguous_array(
            self.coll_ion_cooling_prob_individual
        )

    @property
    def fb_cooling_prob_nd(self):
        return self.fb_cooling_prob_array.shape[1]

    @property
    def coll_exc_cooling_prob_nd(self):
        return self.coll_exc_cooling_prob_array.shape[1]

    @property
    def coll_ion_cooling_prob_nd(self):
        return self.coll_ion_cooling_prob_array.shape[1]
