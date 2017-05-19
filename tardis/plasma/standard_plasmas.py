import os
import logging
import h5py
import pandas as pd
import numpy as np
from tardis.io.config_reader import Configuration
from tardis import atomic
from tardis.util import species_string_to_tuple
from tardis.plasma import BasePlasma
from tardis.plasma.properties.property_collections import (basic_inputs,
    basic_properties, lte_excitation_properties, lte_ionization_properties,
    macro_atom_properties, dilute_lte_excitation_properties,
    nebular_ionization_properties, non_nlte_properties,
    nlte_properties, helium_nlte_properties, helium_numerical_nlte_properties,
    helium_lte_properties, detailed_j_blues_properties, detailed_j_blues_inputs)
from tardis.plasma.exceptions import PlasmaConfigError
from tardis.plasma.properties import (LevelBoltzmannFactorNLTE, JBluesBlackBody,
                                      JBluesDiluteBlackBody, JBluesDetailed,
                                      RadiationFieldCorrection,
                                      StimulatedEmissionFactor,
                                      HeliumNumericalNLTE)

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
        atom_data = atomic.AtomData.from_hdf5(atom_data_fname)

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

    property_kwargs['line_interaction_type'] = config.plasma.line_interaction_type
    property_kwargs['radiative_rates_type'] = config.plasma.radiative_rates_type
    property_kwargs['excitation'] = config.plasma.excitation
    property_kwargs['ionization'] = config.plasma.ionization
    property_kwargs['helium_treatment'] = config.plasma.helium_treatment
    property_kwargs['nlte'] = {}
    property_kwargs['nlte']['species'] = config.plasma.nlte.species

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
        if nlte_conf.classical_nebular and not nlte_conf.coronal_approximation:
            plasma_modules.append(LevelBoltzmannFactorNLTE)
            property_kwargs[LevelBoltzmannFactorNLTE] = {
                'classical_nebular': True}
        elif nlte_conf.coronal_approximation and not nlte_conf.classical_nebular:
            plasma_modules.append(LevelBoltzmannFactorNLTE)
            property_kwargs[LevelBoltzmannFactorNLTE] = {
                'coronal_approximation': True}
        elif nlte_conf.coronal_approximation and nlte_conf.classical_nebular:
            raise PlasmaConfigError('Both coronal approximation and '
                                    'classical nebular specified in the '
                                    'config.')
        property_kwargs[StimulatedEmissionFactor] = dict(
            nlte_species=nlte_species)
    else:
        plasma_modules += non_nlte_properties

    if config.plasma.line_interaction_type in ('downbranch', 'macroatom'):
        plasma_modules += macro_atom_properties

    if config.plasma.helium_treatment == 'recomb-nlte':
        plasma_modules += helium_nlte_properties
        if 'delta_treatment' in config.plasma:
            property_kwargs[RadiationFieldCorrection] = dict(
                delta_treatment=config.plasma.delta_treatment)
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
    kwargs['helium_treatment'] = config.plasma.helium_treatment

    plasma = BasePlasma(plasma_properties=plasma_modules,
                        property_kwargs=property_kwargs, **kwargs)

    return plasma

def from_plasma_hdf(path, h5_file, file_path, model, atom_data_fname=None):
    """
    This function returns a Plasma object 
    from given HDF5 File.

    Parameters
    ----------
    path : 'str'
        Path to transverse in hdf file
    h5_file : 'h5py.File'
        Given HDF5 file
    file_path : 'str'
        Path of Simulation generated HDF file 

    Returns
    -------
    model : `~Plasma`
    """

    if not h5_file:
        raise ValueError("h5_file Parameter can`t be None")
    if atom_data_fname is None:
        atom_data_fname = 'kurucz_cd23_chianti_H_He.h5'

    plasma = {}
    plasma_path = path + '/plasma'
    plasma_dict = {}
    plasma_dict['plasma'] = {}
    plasma_keys = ['property_kwargs', 'metadata']
    atom_data = atomic.AtomData.from_hdf5(atom_data_fname)

    with pd.HDFStore(file_path, 'r') as data:
        for key in h5_file[plasma_path].keys():
            if key in plasma_keys:
                plasma[key] = {}
                buff_path = plasma_path + '/' + key + '/'
                plasma[key] = data[buff_path]

    if plasma['metadata']['atom_data_uuid'] != atom_data.uuid1:
        raise ValueError('Wrong Atom Data passed as parameter')

    property_kwargs = {}
    property_kwargs['plasma'] = {}
    property_kwargs['plasma'] = dict(
        zip(plasma['property_kwargs'].index.format(), plasma['property_kwargs']))
    config = Configuration.from_config_dict(
        property_kwargs, validate=False, plasma_only=True)
    return assemble_plasma(config, model, atom_data=atom_data)