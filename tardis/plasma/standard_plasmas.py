import logging

import pandas as pd
import numpy as np

from tardis.plasma import BasePlasma
from tardis.plasma.properties.property_collections import (basic_inputs,
    basic_properties, lte_excitation_properties, lte_ionization_properties,
    macro_atom_properties, dilute_lte_excitation_properties,
    nebular_ionization_properties, non_nlte_properties,
    nlte_properties, helium_nlte_properties, helium_numerical_nlte_properties)
from tardis.plasma.exceptions import PlasmaConfigError
from tardis.plasma.properties import LevelBoltzmannFactorNLTE

logger = logging.getLogger(__name__)

class LTEPlasma(BasePlasma):

    def __init__(self, t_rad, abundance, density, time_explosion, atomic_data,
        j_blues, link_t_rad_t_electron=0.9, delta_treatment=None):
        plasma_modules = basic_inputs + basic_properties + \
            lte_excitation_properties + lte_ionization_properties + \
            non_nlte_properties

        super(LTEPlasma, self).__init__(plasma_properties=plasma_modules,
            t_rad=t_rad, abundance=abundance, atomic_data=atomic_data,
            density=density, time_explosion=time_explosion, j_blues=j_blues,
	        w=None, link_t_rad_t_electron=link_t_rad_t_electron,
            delta_input=delta_treatment, nlte_species=None)

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

    def initial_w(self, number_densities):
        return np.ones(len(number_densities.columns)) * 0.5

    def update_radiationfield(self, t_rad, ws, j_blues, nlte_config,
        t_electrons=None, n_e_convergence_threshold=0.05,
        initialize_nlte=False):
        if nlte_config is not None and nlte_config.species:
            self.store_previous_properties()
        self.update(t_rad=t_rad, w=ws, j_blues=j_blues)

    def __init__(self, number_densities, atomic_data, time_explosion,
        t_rad=None, delta_treatment=None, nlte_config=None,
        ionization_mode='lte', excitation_mode='lte',
        line_interaction_type='scatter', link_t_rad_t_electron=0.9,
        helium_treatment='none', heating_rate_data_file=None,
        v_inner=None, v_outer=None):

        plasma_modules = basic_inputs + basic_properties

        if excitation_mode == 'lte':
            plasma_modules += lte_excitation_properties
        elif excitation_mode == 'dilute-lte':
            plasma_modules += dilute_lte_excitation_properties
        else:
            raise NotImplementedError('Sorry {0} not implemented yet.'.format(
            excitation_mode))

        if ionization_mode == 'lte':
            plasma_modules += lte_ionization_properties
        elif ionization_mode == 'nebular':
            plasma_modules += nebular_ionization_properties
        else:
            raise NotImplementedError('Sorry ' + ionization_mode +
                ' not implemented yet.')

        if nlte_config is not None and nlte_config.species:
            plasma_modules += nlte_properties
            if nlte_config.classical_nebular==True and \
                nlte_config.coronal_approximation==False:
                LevelBoltzmannFactorNLTE(self, classical_nebular=True)
            elif nlte_config.coronal_approximation==True and \
                nlte_config.classical_nebular==False:
                LevelBoltzmannFactorNLTE(self, coronal_approximation=True)
            elif nlte_config.coronal_approximation==True and \
                nlte_config.classical_nebular==True:
                raise PlasmaConfigError('Both coronal approximation and '
                                        'classical nebular specified in the '
                                        'config.')
        else:
            plasma_modules += non_nlte_properties

        if line_interaction_type in ('downbranch', 'macroatom'):
            plasma_modules += macro_atom_properties

        if t_rad is None:
            t_rad = self.initial_t_rad(number_densities)

        w = self.initial_w(number_densities)

        abundance, density = self.from_number_densities(number_densities,
            atomic_data)

        if nlte_config is not None and nlte_config.species:
            self.nlte_species = nlte_config.species
        else:
            self.nlte_species = None

        if helium_treatment=='recomb-nlte':
            plasma_modules += helium_nlte_properties

        if helium_treatment=='numerical-nlte':
            plasma_modules += helium_numerical_nlte_properties
            if heating_rate_data_file is None:
                raise PlasmaConfigError('Heating rate data file not specified')
            else:
                self.heating_rate_data_file = heating_rate_data_file
                self.v_inner = v_inner
                self.v_outer = v_outer

        self.delta_treatment = delta_treatment

        super(LegacyPlasmaArray, self).__init__(
            plasma_properties=plasma_modules, t_rad=t_rad,
            abundance=abundance, density=density,
            atomic_data=atomic_data, time_explosion=time_explosion,
            j_blues=None, w=w, link_t_rad_t_electron=link_t_rad_t_electron)
