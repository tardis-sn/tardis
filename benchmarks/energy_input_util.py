"""
Basic TARDIS Benchmark.
"""
import numpy as np
from asv_runner.benchmarks.mark import parameterize, skip_benchmark

import tardis.energy_input.util as util
from benchmarks.benchmark_base import BenchmarkBase
from tardis.energy_input.util import R_ELECTRON_SQUARED, get_perpendicular_vector


# @skip_benchmark
class BenchmarkEnergyInputUtil(BenchmarkBase):
    """
    Class to benchmark the util function.
    """

    def __init__(self):
        pass

    @parameterize({"Spherical to cartesian": [
        {
            "r": 1,
            "theta": 0,
            "phi": 0,
        },
        {
            "r": 1,
            "theta": np.pi,
            "phi": 0,
        },
        {
            "r": 1,
            "theta": np.pi / 2,
            "phi": 0,
        },
        {
            "r": 1,
            "theta": np.pi / 2,
            "phi": np.pi,
        },
        {
            "r": 1,
            "theta": np.pi / 2,
            "phi": np.pi / 2,
        },
    ]})
    def time_spherical_to_cartesian(self, values: dict):
        util.spherical_to_cartesian(values.get('r'), values.get('theta'), values.get('phi'))

    @parameterize({"Kappa calculation": [
        {
            "energy": 511.0,
        },
        {
            "energy": 255.5,
        },
        {
            "energy": 0.0,
        },
        {
            "energy": 511.0e7,
        },
    ]})
    def time_kappa_calculation(self, values: dict):
        util.kappa_calculation(values.get('energy'))

    @parameterize({"Klein Nishina": [
        {
            "energy": 511.0e3,
            "theta_C": 1.0
        },
        {
            "energy": 255.5e3,
            "theta_C": np.pi
        },
        {
            "energy": 0.0,
            "theta_C": 2.0 * np.pi
        },
        {
            "energy": 511.0e10,
            "theta_C": np.pi / 2.0
        },
    ]})
    def time_klein_nishina(self, values: dict):
        util.klein_nishina(values.get('energy'), values.get('theta_C'))

        kappa = util.kappa_calculation(values.get('energy'))

        (
                R_ELECTRON_SQUARED
                / 2
                * (1.0 + kappa * (1.0 - np.cos(values.get('theta_C')))) ** -2.0
                * (
                        1.0
                        + np.cos(values.get('theta_C')) ** 2.0
                        + (kappa ** 2.0 * (1.0 - np.cos(values.get('theta_C'))) ** 2.0)
                        / (1.0 + kappa * (1.0 - np.cos(values.get('theta_C'))))
                )
        )

    def time_get_perpendicular_vector(self):
        input_vector = np.array([0.3, 0.4, 0.5])
        get_perpendicular_vector(input_vector)
