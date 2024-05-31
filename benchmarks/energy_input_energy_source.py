from asv_runner.benchmarks.mark import parameterize
import numpy as np
import pandas as pd
import radioactivedecay as rd

import tardis.energy_input.energy_source as energy_source
from benchmarks.benchmark_base import BenchmarkBase

from tardis.util.base import (
    atomic_number2element_symbol,
)
from tardis.energy_input.util import (
    convert_half_life_to_astropy_units,
    ELECTRON_MASS_ENERGY_KEV,
)


class BenchmarkEnergyInputEnergySource(BenchmarkBase):
    """
    Class to benchmark the energy source
    """

    def time_positronium_continuum(self):
        energy_source.positronium_continuum()

    @parameterize(
        {
            "atomic number": [
                1, 
                24, 
                13, 
                18,
            ], 
            "atomic mass": [
                1.008, 
                51.9961,
                26.9815385,
                39.948,
            ]
        }
    )
    def time_get_isotope_string(atomic_number, atomic_mass):
        energy_source.get_isotope_string(atomic_number, atomic_mass)