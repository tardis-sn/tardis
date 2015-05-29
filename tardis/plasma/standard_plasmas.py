import logging

from tardis.plasma import BasePlasma
from tardis.plasma.properties.property_collections import (basic_inputs,
    lte_processing_properties)

logger = logging.getLogger(__name__)

class LTEPlasma(BasePlasma):

    def __init__(self, t_rad, abundance, density, time_explosion, atomic_data):
        plasma_modules = basic_inputs + lte_processing_properties

        super(LTEPlasma, self).__init__(plasma_modules=plasma_modules,
                                        t_rad=t_rad, abundance=abundance,
                                        atomic_data=atomic_data,
                                        density=density,
                                        time_explosion=time_explosion)

class LegacyPlasma(BasePlasma):

    def __init__(self, number_densities, atom_data, time_explosion,
                 delta_treatment=None, nlte_config=None, ionization_mode='lte',
                 excitation_mode='lte'):
        pass
