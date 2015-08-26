import logging
import pandas as pd
import numpy as np

from tardis.plasma import BasePlasma
from tardis.plasma.properties.property_collections import (basic_inputs,
    basic_properties, lte_excitation_properties, lte_ionization_properties,
    macro_atom_properties, dilute_lte_excitation_properties,
    nebular_ionization_properties, non_nlte_properties,
    nlte_properties, helium_nlte_properties, nlte_ionization_properties,
    non_nlte_ionization_properties)
from tardis.io.util import parse_abundance_dict_to_dataframe

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
            delta_input=delta_treatment, nlte_excitation_species=None,
            nlte_ionization_species=None, previous_beta_sobolevs=None,
            previous_electron_densities=None)

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

    def update_radiationfield(self, t_rad, ws, j_blues, nlte_excitation_config,
        t_electrons=None, n_e_convergence_threshold=0.05,
        initialize_nlte_excitation=False):
        if (nlte_excitation_config is not None and
            nlte_excitation_config.species):
            self.store_previous_properties()
        self.update(t_rad=t_rad, w=ws, j_blues=j_blues)

    def __init__(self, number_densities, atomic_data, time_explosion,
        t_rad=None, delta_treatment=None, nlte_excitation_config=None,
        nlte_ionization_config=None, ionization_mode='lte',
        excitation_mode='lte', line_interaction_type='scatter',
        link_t_rad_t_electron=0.9, helium_treatment='lte'):

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

        if (nlte_excitation_config is not None and
            nlte_excitation_config.species):
            plasma_modules += nlte_properties
            nlte_excitation_species = nlte_excitation_config.species
        else:
            plasma_modules += non_nlte_properties
            nlte_excitation_species = None

        if (nlte_ionization_config is not None and
            nlte_ionization_config.species):
            plasma_modules += nlte_ionization_properties
            nlte_ionization_species = nlte_ionization_config.species
        else:
            plasma_modules += non_nlte_ionization_properties
            nlte_ionization_species = None

        if line_interaction_type in ('downbranch', 'macroatom'):
            plasma_modules += macro_atom_properties

        if t_rad is None:
            t_rad = self.initial_t_rad(number_densities)

        w = self.initial_w(number_densities)

        abundance, density = self.from_number_densities(number_densities,
            atomic_data)

        try:
            initial_beta_sobolevs = np.ones((len(atomic_data.lines),
                len(number_densities.columns)))
        except:
            initial_beta_sobolevs = np.ones((len(atomic_data._lines),
                len(number_densities.columns)))
        initial_electron_densities = number_densities.sum(axis=0)

        self.nlte_excitation_config = nlte_excitation_config

        self.nlte_ionization_config = nlte_ionization_config

        if helium_treatment=='recomb-nlte':
            plasma_modules += helium_nlte_properties

        super(LegacyPlasmaArray, self).__init__(plasma_properties=plasma_modules,
            t_rad=t_rad, abundance=abundance, density=density,
            atomic_data=atomic_data, time_explosion=time_explosion,
            j_blues=None, w=w, link_t_rad_t_electron=link_t_rad_t_electron,
            delta_input=delta_treatment,
            nlte_excitation_species=nlte_excitation_species,
            nlte_ionization_species=nlte_ionization_species,
            previous_electron_densities=initial_electron_densities,
            previous_beta_sobolevs=initial_beta_sobolevs)
