import logging

import pandas as pd
import numpy as np

from tardis.plasma import BasePlasma
from tardis.plasma.properties.property_collections import (basic_inputs,
    basic_properties, lte_excitation_properties, lte_ionization_properties,
    macro_atom_properties, dilute_lte_excitation_properties,
    nebular_ionization_properties, non_nlte_properties,
    nlte_properties, helium_nlte_properties, helium_numerical_nlte_properties,
    helium_lte_properties)
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


def assemble_plasma(config, model):
    """
    Create a BasePlasma instance from a Configuration object
    and a Radial1DModel.

    Parameters
    ----------
    config: ~io.config_reader.Configuration
    model: ~model.Radial1DModel

    Returns
    -------
    : ~plasma.BasePlasma

    """
    kwargs = dict(t_rad=model.t_radiative, abundance=model.abundance,
                  density=model.density, atomic_data=config.atom_data,
                  time_explosion=model.time_explosion, j_blues=None,
                  w=model.dilution_factor, link_t_rad_t_electron=0.9)

    plasma_modules = basic_inputs + basic_properties

    if config.plasma.excitation == 'lte':
        plasma_modules += lte_excitation_properties
    elif config.plasma.excitation == 'dilute-lte':
        plasma_modules += dilute_lte_excitation_properties

    if config.plasma.ionization == 'lte':
        plasma_modules += lte_ionization_properties
    elif config.plasma.ionization == 'nebular':
        plasma_modules += nebular_ionization_properties

    nlte_conf = config.plasma.nlte
    if nlte_conf.species:
        plasma_modules += nlte_properties
        kwargs['nlte_species'] = nlte_conf.species
        if nlte_conf.classical_nebular and not nlte_conf.coronal_approximation:
            plasma_modules.append(LevelBoltzmannFactorNLTE)
            kwargs['classical_nebular'] = True
        elif nlte_conf.coronal_approximation and not nlte_conf.classical_nebular:
            plasma_modules.append(LevelBoltzmannFactorNLTE)
            kwargs['coronal_approximation'] = True
        elif nlte_conf.coronal_approximation and nlte_conf.classical_nebular:
            raise PlasmaConfigError('Both coronal approximation and '
                                    'classical nebular specified in the '
                                    'config.')
    else:
        plasma_modules += non_nlte_properties

    if config.plasma.line_interaction_type in ('downbranch', 'macroatom'):
        plasma_modules += macro_atom_properties

    if config.plasma.helium_treatment == 'recomb-nlte':
        plasma_modules += helium_nlte_properties
    elif config.plasma.helium_treatment == 'numerical-nlte':
        plasma_modules += helium_numerical_nlte_properties
        # TODO: See issue #633
        if config.plasma.heating_rate_data_file in ['none', None]:
            raise PlasmaConfigError('Heating rate data file not specified')
        else:
            kwargs['heating_rate_data_file'] = \
                config.plasma.heating_rate_data_file
            kwargs['v_inner'] = model.v_inner
            kwargs['v_outer'] = model.v_outer
    else:
        plasma_modules += helium_lte_properties
    kwargs['helium_treatment'] = config.plasma.helium_treatment

    if 'delta_treatment' in config.plasma:
        kwargs['delta_treatment'] = config.plasma.delta_treatment

    plasma = BasePlasma(plasma_properties=plasma_modules, **kwargs)

    return plasma
