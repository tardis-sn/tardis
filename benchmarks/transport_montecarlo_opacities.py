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
            "Parameters": [
                {
                    "Ejecta_density": 0.01,
                    "Energy": 511.0,
                    "Iron_group_fraction": 0.5,
                },
                {
                    "Ejecta_density": 0.01,
                    "Energy": 5110000000.0,
                    "Iron_group_fraction": 0.0,
                },
                {
                    "Ejecta_density": 0.01,
                    "Energy": 255.5,
                    "Iron_group_fraction": 0.0,
                },
                {
                    "Ejecta_density": 0.01,
                    "Energy": 255.5,
                    "Iron_group_fraction": 0.5,
                },
                {
                    "Ejecta_density": 100000.0,
                    "Energy": 255.5,
                    "Iron_group_fraction": 1.0,
                }
            ]
        }
    )
    def time_photoabsorption_opacity_calculation(self, parameters):
        calculate_opacity.photoabsorption_opacity_calculation(
            parameters["Energy"],
            parameters["Ejecta_density"],
            parameters["Iron_group_fraction"],
        )

    @parameterize(
        {
            "Parameters": [
                {
                    "Ejecta_density": 0.01,
                    "Energy": 511.0,
                    "Iron_group_fraction": 0.5,
                },
                {
                    "Ejecta_density": 0.01,
                    "Energy": 5110000000.0,
                    "Iron_group_fraction": 0.0,
                },
                {
                    "Ejecta_density": 0.01,
                    "Energy": 255.5,
                    "Iron_group_fraction": 0.0,
                },
                {
                    "Ejecta_density": 0.01,
                    "Energy": 255.5,
                    "Iron_group_fraction": 0.5,
                },
                {
                    "Ejecta_density": 100000.0,
                    "Energy": 255.5,
                    "Iron_group_fraction": 1.0,
                }
            ]
        }
    )
    def time_pair_creation_opacity_calculation(self, parameters):
        calculate_opacity.pair_creation_opacity_calculation(
            parameters["Energy"],
            parameters["Ejecta_density"],
            parameters["Iron_group_fraction"],
        )
