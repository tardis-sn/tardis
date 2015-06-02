import logging
import pandas as pd
import numpy as np

from tardis.plasma import BasePlasma
from tardis.plasma.properties.property_collections import (basic_inputs,
    basic_properties, lte_excitation_properties, lte_ionization_properties)

logger = logging.getLogger(__name__)

class LTEPlasma(BasePlasma):

    def __init__(self, t_rad, abundance, density, time_explosion, atomic_data):
        plasma_modules = basic_inputs + basic_properties + \
            lte_excitation_properties + lte_ionization_properties

        super(LTEPlasma, self).__init__(plasma_modules=plasma_modules,
            t_rad=t_rad, abundance=abundance, atomic_data=atomic_data,
            density=density, time_explosion=time_explosion)

class LegacyPlasmaArray(BasePlasma):

    def from_number_densities(self, number_densities, atomic_data):
        atomic_mass = atomic_data.atom_data.ix[number_densities.index].mass
        elemental_density = number_densities.mul(atomic_mass,
            axis='index')
        density = elemental_density.sum()
        abundance = pd.DataFrame(elemental_density/density,
            index=number_densities.index, columns=number_densities.columns,
            dtype=np.float64)
        return abundance, density

    def initial_t_rad(self, number_densities):
        return np.ones(len(number_densities.columns)) * 10000

    def update_radiationfield(self, t_rad, ws, j_blues,
        t_electrons=None, n_e_convergence_threshold=0.05,
        initialize_nlte=False):
        self.update(t_rad=t_rad)

    def __init__(self, number_densities, atomic_data, time_explosion,
        t_rad=None, delta_treatment=None, nlte_config=None,
        ionization_mode='lte', excitation_mode='lte'):

        plasma_modules = basic_inputs + basic_properties
        if excitation_mode == 'lte':
            plasma_modules += lte_excitation_properties
        elif excitation_mode == 'dilute-lte':
            raise NotImplementedError('Sorry ' + excitation_mode +
                ' not implemented yet.')
        if ionization_mode == 'lte':
            plasma_modules += lte_ionization_properties
        elif ionization_mode == 'nebular':
            raise NotImplementedError('Sorry ' + ionization_mode +
                ' not implemented yet.')
        if nlte_config.species:
            raise NotImplementedError('Sorry, NLTE treatment not implemented \
                yet.')
        if delta_treatment is not None:
            raise NotImplementedError('Sorry, delta treatment not implemented \
                yet')

        if t_rad is None:
            t_rad = self.initial_t_rad(number_densities)

        abundance, density = self.from_number_densities(number_densities,
            atomic_data)

        super(LegacyPlasmaArray, self).__init__(plasma_modules=plasma_modules,
            t_rad=t_rad, abundance=abundance, density=density,
            atomic_data=atomic_data, time_explosion=time_explosion)
