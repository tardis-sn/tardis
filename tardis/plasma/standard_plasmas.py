import os
import logging

import pandas as pd
import numpy as np

from tardis.io import atomic
from tardis.util.base import species_string_to_tuple
from tardis.plasma import BasePlasma
from tardis.plasma.properties.property_collections import (basic_inputs,
    basic_properties, lte_excitation_properties, lte_ionization_properties,
    macro_atom_properties, dilute_lte_excitation_properties,
    nebular_ionization_properties, non_nlte_properties,
    nlte_properties, helium_nlte_properties, helium_numerical_nlte_properties,
    helium_lte_properties, detailed_j_blues_properties, detailed_j_blues_inputs)
from tardis.plasma.exceptions import PlasmaConfigError

from tardis.plasma.properties import (
        LevelBoltzmannFactorNLTE,
        JBluesBlackBody,
        JBluesDiluteBlackBody,
        JBluesDetailed,
        RadiationFieldCorrection,
        StimulatedEmissionFactor,
        HeliumNumericalNLTE,
        IonNumberDensity)

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


def assemble_plasma(config, model, atom_data=None):
    """
    Create a BasePlasma instance from a Configuration object
    and a Radial1DModel.

    Parameters
    ----------
    config: ~io.config_reader.Configuration
    model: ~model.Radial1DModel
    atom_data: ~atomic.AtomData
        If None, an attempt will be made to read the atomic data
        from config.

    Returns
    -------
    : ~plasma.BasePlasma

    """
    # Convert the nlte species list to a proper format.
    nlte_species = [species_string_to_tuple(s) for s in
                    config.plasma.nlte.species]

    if atom_data is None:
        if 'atom_data' in config:
            if os.path.isabs(config.atom_data):
                atom_data_fname = config.atom_data
            else:
                atom_data_fname = os.path.join(config.config_dirname,
                                               config.atom_data)
        else:
            raise ValueError('No atom_data option found in the configuration.')

        logger.info('Reading Atomic Data from %s', atom_data_fname)

        try:
            atom_data = atomic.AtomData.from_hdf(atom_data_fname)
        except TypeError as e:
            print(e, 'Error might be from the use of an old-format of the atomic database, \n'
                'please see https://github.com/tardis-sn/tardis-refdata/tree/master/atom_data'
                ',for the most recent version.')
            raise

    atom_data.prepare_atom_data(
        model.abundance.index,
        line_interaction_type=config.plasma.line_interaction_type,
        nlte_species=nlte_species)

    kwargs = dict(t_rad=model.t_radiative, abundance=model.abundance,
                  density=model.density, atomic_data=atom_data,
                  time_explosion=model.time_explosion,
                  w=model.dilution_factor, link_t_rad_t_electron=0.9)

    plasma_modules = basic_inputs + basic_properties
    property_kwargs = {}
    if config.plasma.radiative_rates_type == 'blackbody':
        plasma_modules.append(JBluesBlackBody)
    elif config.plasma.radiative_rates_type == 'dilute-blackbody':
        plasma_modules.append(JBluesDiluteBlackBody)
    elif config.plasma.radiative_rates_type == 'detailed':
        plasma_modules += detailed_j_blues_properties + detailed_j_blues_inputs
        kwargs.update(r_inner=model.r_inner,
                      t_inner=model.t_inner,
                      volume=model.volume,
                      j_blue_estimator=None)
        property_kwargs[JBluesDetailed] = {'w_epsilon': config.plasma.w_epsilon}
    else:
        raise ValueError('radiative_rates_type type unknown - %s',
                         config.plasma.radiative_rates_type)

    if config.plasma.excitation == 'lte':
        plasma_modules += lte_excitation_properties
    elif config.plasma.excitation == 'dilute-lte':
        plasma_modules += dilute_lte_excitation_properties

    if config.plasma.ionization == 'lte':
        plasma_modules += lte_ionization_properties
    elif config.plasma.ionization == 'nebular':
        plasma_modules += nebular_ionization_properties

    if nlte_species:
        plasma_modules += nlte_properties
        nlte_conf = config.plasma.nlte
        plasma_modules.append(
                LevelBoltzmannFactorNLTE.from_config(nlte_conf)
                )
        property_kwargs[StimulatedEmissionFactor] = dict(
            nlte_species=nlte_species)
    else:
        plasma_modules += non_nlte_properties

    if config.plasma.line_interaction_type in ('downbranch', 'macroatom'):
        plasma_modules += macro_atom_properties

    if 'delta_treatment' in config.plasma:
        property_kwargs[RadiationFieldCorrection] = dict(
            delta_treatment=config.plasma.delta_treatment)

    if config.plasma.helium_treatment == 'recomb-nlte':
        plasma_modules += helium_nlte_properties
    elif config.plasma.helium_treatment == 'numerical-nlte':
        plasma_modules += helium_numerical_nlte_properties
        # TODO: See issue #633
        if config.plasma.heating_rate_data_file in ['none', None]:
            raise PlasmaConfigError('Heating rate data file not specified')
        else:
            property_kwargs[HeliumNumericalNLTE] = dict(
                heating_rate_data_file=config.plasma.heating_rate_data_file)
    else:
        plasma_modules += helium_lte_properties

    if model._electron_densities:
        electron_densities = pd.Series(model._electron_densities.cgs.value)
        if config.plasma.helium_treatment == 'numerical-nlte':
            property_kwargs[IonNumberDensityHeNLTE] = dict(
                electron_densities=electron_densities)
        else:
            property_kwargs[IonNumberDensity] = dict(
                electron_densities=electron_densities)

    kwargs['helium_treatment'] = config.plasma.helium_treatment

    plasma = BasePlasma(plasma_properties=plasma_modules,
                        property_kwargs=property_kwargs, **kwargs)

    return plasma
