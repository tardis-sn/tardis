"""
Basic TARDIS Benchmark.
"""

import tardis.opacities.opacities as calculate_opacity
from benchmarks.benchmark_base import BenchmarkBase
from tardis.opacities.opacities import compton_opacity_calculation


class BenchmarkMontecarloMontecarloNumbaOpacities(BenchmarkBase):
    """
    Class to benchmark the numba opacities function.
    """

    def time_compton_opacity_calculation(self):
        energy = 511.0
        electron_number_density = 1e15
        calculate_opacity.compton_opacity_calculation(
            energy, electron_number_density
        )

    def time_photoabsorption_opacity_calculation(self):
        energy = 255.5
        ejecta_density = 100000.0
        iron_group_fraction = 0.5
        calculate_opacity.photoabsorption_opacity_calculation(
            energy, ejecta_density, iron_group_fraction
        )

    def time_pair_creation_opacity_calculation(self):
        energy = 255.9
        ejecta_density = 100000.009
        iron_group_fraction = 0.5
        calculate_opacity.pair_creation_opacity_calculation(
            energy, ejecta_density, iron_group_fraction
        )
