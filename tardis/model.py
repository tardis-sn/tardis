#File contains model class

import numpy as np
import os
import plasma
import logging
import constants
import texttable

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
        time_0 = 0.000231481 * constants.days2seconds
        densities = 10 ** densities * (current_time / time_0) ** -3
        #TODO implement abundances here.
        return cls(velocities[112::5], densities[112::5], current_time)

    @classmethod
    def from_lucy99(cls, v_inner, abundances, current_time, no_of_shells=20, v_outer=30000 * 1e5,
                    abundance_mode='uniform', density_coefficient=3e29):
        """
        :param cls:
        :param v_inner: in km/s
        :param current_time: in days
        :param shell_velocity: in km/s
        :param no_of_shells: number of shells
        :param density_coefficient: see lucy 99
        :param abundance_mode: currently supported uniform, scaled
        :return:
        """

        time_0 = 0.000231481 * constants.days2seconds

        #TODO are the shells distributed right? currently just linear.
        velocities = np.linspace(v_inner, v_outer, no_of_shells, endpoint=True)
        densities = density_coefficient * (velocities * 1e-5) ** -7
        densities *= (current_time / time_0) ** -3

        if abundance_mode == 'uniform':
            model_abundances = np.empty((no_of_shells, abundances.size))
            for i in xrange(no_of_shells):
                model_abundances[i] = abundances.copy()

        elif abundance_mode == 'scaled':
            raise NotImplementedError()

        return cls(velocities, densities, model_abundances, current_time)

    @classmethod
    def vstruct_exponential(cls, roh0, v_inner, abundances, current_time, no_of_shells=20, v_outer=30000 * 1e5 ):
    #Legacy
        densities = roh0 * (current_time / time_0) ** -3

    def __init__(self, velocities, densities, abundances, current_time, ws=None):
        """
        :param velocities: in cm/s
        :param densities: in g/cm^3 (I think)
        :param current_time: in seconds
        """
        self.v_inner = velocities[:-1]
        self.v_outer = velocities[1:]

        self.current_time = current_time

        self.r_inner = self.v_inner * current_time
        self.r_outer = self.v_outer * current_time
        self.r_middle = 0.5 * (self.r_inner + self.r_outer)
        if ws is None:
            self.ws = 0.5 * (1 - np.sqrt(1 - self.r_inner[0] ** 2 / self.r_middle ** 2))
        else:
            self.ws = np.array([(0.5 * (1 - np.sqrt(1 - self.r_inner[0] ** 2 / self.r_middle[i] ** 2))) if w < 0\
                                else w for i, w in enumerate(ws)])

        self.densities_middle = densities[1:]
        print "sdsadsad", self.densities_middle
        self.abundances = abundances


    def set_atomic_model(self, atomic_model):
        self.atomic_model = atomic_model


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

    def initialize_plasmas(self, t_rads):
        #self._initialize_temperature(t_rad)
        self.t_rads = t_rads
        self.electron_densities = np.empty(self.r_inner.size, dtype=np.float64)
        self.plasmas = []
        for i, (current_abundance, current_t, current_w, current_density) in enumerate(
            zip(self.abundances, self.t_rads, self.ws, self.densities_middle)):
            current_plasma = plasma.NebularPlasma.from_model(current_abundance, current_density, self.atomic_model)
            current_plasma.update_radiationfield(t_rad=current_t, w=current_w)
            self.plasmas.append(current_plasma)
            self.electron_densities[i] = current_plasma.electron_density

    def print_model_table(self):
        temp_table = texttable.Texttable()
        header = ('v_inner', 'v_outer', 'r_inner', 'r_outer', 'density', 't_rad', 'w')
        temp_table.add_row(header)
        header = ('[km/s]', '[km/s]', '[km]', '[km]', '[g/cm^3]', '[K]', '[1]')
        temp_table.add_row(header)
        temp_table.set_deco(temp_table.HEADER | temp_table.VLINES)
        #TODO strangely enough '%.2e' doesn't seem to work for densities. needs to set a dot.
        for v_inner, v_outer, r_inner, r_outer, density, t_rad, w in zip(self.v_inner, self.v_outer, self.r_inner,
            self.r_outer, self.densities_middle,
            self.t_rads, self.ws):
            temp_table.add_row(('%.2f' % (v_inner / 1e5,), '%.2f' % (v_outer / 1e5,), '%.2e' % (r_inner / 1e5,),
                                '%.2e' % (r_outer / 1e5,), '%.2e .' % density, '%.2f' % t_rad, '%.2f' % w))
        return temp_table.draw()

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