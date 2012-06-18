#File contains model class

import numpy as np
import os
import plasma
import logging

logger = logging.getLogger(__name__)

miss_distance = 1e99

#TODO use pkgutil to read data files
w7model_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data', 'w7model.dat'))
w7abundances_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data', 'w7abundances.dat'))

class Model(object):
    pass


class MultiZoneRadial(Model):
    @classmethod
    def from_w7(cls, current_time):
        velocities, densities = np.loadtxt(w7model_file, usecols=(1, 2), unpack=1)
        velocities *= 1e5
        densities = 10 ** densities
        time_0 = 0.000231481 * 86400
        return cls(velocities[112::5], densities[112::5], current_time, time_0)

    @classmethod
    def from_lucy99(cls, v_inner, current_time, shell_velocity=386, shells=20, density_coefficient=3e29):
        """

        :param cls:
        :param v_inner: in km/s
        :param current_time: in days
        :param shell_velocity: in km/s
        :param shells: number of shells
        :param density_coefficient: see lucy 99
        :return:
        """
        time_0 = 0.000231481 * 86400
        velocities = np.arange(shells, dtype=np.float64) * shell_velocity + v_inner
        densities = density_coefficient * velocities ** -7
        velocities *= 1e5
        return cls(velocities, densities, current_time, time_0)

    def __init__(self, velocities, densities, current_time, time_0):
        #constructing zones
        self.v_inner = velocities[:-1]
        self.v_outer = velocities[1:]

        self.time_0 = time_0
        self.current_time = current_time

        self.r_inner = self.v_inner * current_time
        self.r_outer = self.v_outer * current_time
        self.r_middle = 0.5 * (self.r_inner + self.r_outer)
        self.ws = 0.5 * (1 - np.sqrt(1 - self.r_inner[0] ** 2 / self.r_middle ** 2))
        self.densities_middle = (current_time / time_0) ** -3 * densities[1:]

        if logger.getEffectiveLevel() == logging.DEBUG:
            model_string = "New Model:\n%s\n" % ('-' * 80,)
            model_string += "% 10s% 10s% 10s% 10s% 10s\n" %\
                            ('v_inner', 'v_outer', 'r_inner', 'r_outer', 'density')

            model_string += "% 10s% 10s% 10s% 10s% 10s\n" %\
                            ('[km/s]', '[km/s]', '[km]', '[km]', '[g/cm^3]')
            for line in zip(self.v_inner,
                self.v_outer,
                self.r_inner,
                self.r_outer,
                self.densities_middle):
                model_string += "% 10.2f% 10.2f% 10.2e% 10.2e% 10.2e\n"
                model_string = model_string % (line[0] / 1e5, line[1] / 1e5, line[2] / 1e5, line[3] / 1e5, line[4])

            logger.debug(model_string)


    def set_atomic_model(self, atomic_model):
        self.atomic_model = atomic_model

    def read_abundances_uniform(self, abundances):
        self.abundances = np.empty((self.r_inner.size, abundances.size))
        for i in range(self.r_inner.size):
            self.abundances[i] = abundances.copy()

    def read_w7_abundances(self):
        w7_abundances = np.loadtxt(w7abundances_file, usecols=np.arange(1, 31))
        w7ni56_0 = np.loadtxt(w7model_file, usecols=(4,))
        ni56, co56, fe56 = set_ni_chain(self.current_time, w7ni56_0)
        w7_abundances[:, 27] += - w7ni56_0 + ni56
        w7_abundances[:, 26] += co56
        w7_abundances[:, 25] += fe56
        self.abundances = w7_abundances[90::6]
        #NORMALIZING
        self.abundances.sum(axis=1).reshape(self.abundances.shape[0], 1)

    def _initialize_temperature(self, t_rad):
        self.t_rads = np.ones(self.r_inner.size) * t_rad

    def initialize_plasmas(self, t_rad):
        self._initialize_temperature(t_rad)
        self.electron_densities = np.empty(self.r_inner.size, dtype=np.float64)
        self.plasmas = []
        for i, (current_abundance, current_t, current_w, current_density) in enumerate(
            zip(self.abundances, self.t_rads, self.ws, self.densities_middle)):
            current_plasma = plasma.NebularPlasma.from_model(current_abundance, current_density, self.atomic_model,
                abundance_mode='array')
            current_plasma.update_radiationfield(t_rad=current_t, w=current_w)
            self.plasmas.append(
                current_plasma
            )
            self.electron_densities[i] = current_plasma.electron_density


    def calculate_tau_sobolevs(self):
        tau_sobolevs = np.empty((self.r_inner.size, len(self.atomic_model.line_list)))
        for i, plasma in enumerate(self.plasmas):
            tau_sobolevs[i] = plasma.calculate_tau_sobolev(self.current_time)
        return tau_sobolevs

    def calculate_transition_probabilities(self, tau_sobolevs):
        transition_probabilities = np.empty((self.r_inner.size, self.atomic_model.macro_atom.target_level_total.size))
        for i, (t_rad, w, tau_sobolev) in enumerate(zip(self.t_rads, self.ws, tau_sobolevs)):
            transition_probabilities[i] = self.atomic_model.macro_atom.calculate_transition_probabilities(tau_sobolev,
                t_rad, w)
        return transition_probabilities

    def update_model(self, t_rads, ws):
        self.t_rads = t_rads
        self.ws = ws
        for i, (w, t_rad, current_plasma) in enumerate(zip(self.ws, self.t_rads, self.plasmas)):
            current_plasma.update_radiationfield(t_rad, w)
            self.electron_densities[i] = current_plasma.electron_density


def set_ni_chain(time_exp, ni56):
    t1 = 6.1
    t2 = 77.27
    ni_abundance_0 = ni56
    ni_abundance_t = ni_abundance_0 * (pow(2, -(time_exp / t1)))
    co_abundance_t = (t2 / (t1 - t2)) * ni_abundance_0 * (pow(2, -time_exp / t1) - pow(2, -time_exp / t2))
    fe_abundance_t = ni_abundance_0 - co_abundance_t - ni_abundance_t
    return ni_abundance_t, co_abundance_t, fe_abundance_t