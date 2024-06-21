"""
Basic TARDIS Benchmark.
"""

from asv_runner.benchmarks.mark import parameterize

import tardis.opacities.opacities as calculate_opacity
from benchmarks.benchmark_base import BenchmarkBase


class BenchmarkMontecarloMontecarloNumbaOpacities(BenchmarkBase):
    """
    Class to benchmark the numba opacities function.
    """

    @parameterize(
        {
            "Electron number density": [
                1.0e11,
                1e15,
                1e5,
            ],
            "Energy": [
                511.0,
                255.5,
                511.0e7,
            ],
        }
    )
    def time_compton_opacity_calculation(self, electron_number_density, energy):
        calculate_opacity.compton_opacity_calculation(
            energy, electron_number_density
        )

    @parameterize(
        {
            "Ejecta density": [
                1.0,
                1e-2,
                1e-2,
                1e5,
            ],
            "Energy": [
                511.0,
                255.5,
                255.5,
                511.0e7,
            ],
            "Iron group fraction": [
                0.0,
                0.5,
                0.25,
                1.0,
            ],
        }
    )
    def time_photoabsorption_opacity_calculation(
        self, ejecta_density, energy, iron_group_fraction
    ):
        calculate_opacity.photoabsorption_opacity_calculation(
            energy, ejecta_density, iron_group_fraction
        )

    @parameterize(
        {
            "Ejecta density": [
                1.0,
                1e-2,
                1e-2,
                1e5,
            ],
            "Energy": [
                511.0,
                1500,
                1200,
                511.0e7,
            ],
            "Iron group fraction": [
                0.0,
                0.5,
                0.25,
                1.0,
            ],
        }
    )
    def time_pair_creation_opacity_calculation(
        self, ejecta_density, energy, iron_group_fraction
    ):
        calculate_opacity.pair_creation_opacity_calculation(
            energy, ejecta_density, iron_group_fraction
        )
