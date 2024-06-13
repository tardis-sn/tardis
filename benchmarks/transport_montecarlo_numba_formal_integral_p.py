"""
Basic TARDIS Benchmark.
"""

import numpy as np
from asv_runner.benchmarks.mark import parameterize, skip_benchmark

import tardis.transport.montecarlo.formal_integral as formal_integral
from benchmarks.benchmark_base import BenchmarkBase
from tardis import constants as c
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.util.base import intensity_black_body


class BenchmarkMontecarloMontecarloNumbaNumbaFormalIntegral(BenchmarkBase):
    """
    Class to benchmark the numba formal integral function.
    """

    @parameterize(
        {
            "nu": [1e14, 0, 1],
            "temperature": [1e4, 1, 1],
        }
    )
    def time_intensity_black_body(self, nu, temperature):
        func = formal_integral.intensity_black_body
        actual = func(nu, temperature)
        print(actual, type(actual))
        intensity_black_body(nu, temperature)

    @parameterize({"N": (1e2, 1e3, 1e4, 1e5)})
    def time_trapezoid_integration(self, n):
        func = formal_integral.trapezoid_integration
        h = 1.0
        n = int(n)
        data = np.random.random(n)

        func(data, h)
        np.trapz(data)

    @staticmethod
    def calculate_z(r, p):
        return np.sqrt(r * r - p * p)

    TESTDATA = [
        np.linspace(1, 2, 3, dtype=np.float64),
        np.linspace(0, 1, 3),
        # np.linspace(1, 2, 10, dtype=np.float64),
    ]

    def formal_integral_geometry(self, r):
        # NOTE: PyTest is generating a full matrix with all the permutations.
        #       For the `time_calculate_z` function with values: [0.0, 0.5, 1.0]
        #       -  p=0.0, formal_integral_geometry0-0.0, param["r"]: [1.  1.5 2. ]
        #       -  p=0.5, formal_integral_geometry0-0.5, param["r"]: [1.  1.5 2. ]
        #       -  p=1.0, formal_integral_geometry0-1.0, param["r"]: [1.  1.5 2. ]
        #       -  p=0.0, formal_integral_geometry1-0.0, param["r"]: [0.  0.5 1. ]
        #       -  p=1.0, formal_integral_geometry1-1.0, param["r"]: [0.  0.5 1. ]
        #       Same for `test_populate_z_photosphere` function
        #       And for `test_populate_z_shells` function
        #       -  p=1e-05, formal_integral_geometry0-1e-05, param["r"]: [1.  1.5 2. ]
        #       -  p=0.5,   formal_integral_geometry0-0.5,   param["r"]: [1.  1.5 2. ]
        #       -  p=0.99,  formal_integral_geometry0-0.99,  param["r"]: [1.  1.5 2. ]
        #       -  p=1,     formal_integral_geometry0-1,     param["r"]: [1.  1.5 2. ]
        #       -  p=1e-05, formal_integral_geometry1-1e-05, param["r"]: [0.  0.5 1. ]
        #       -  p=0.5,   formal_integral_geometry1-0.5,   param["r"]: [0.  0.5 1. ]
        #       -  p=0.99,  formal_integral_geometry1-0.99,  param["r"]: [0.  0.5 1. ]
        #       -  p=1,     formal_integral_geometry1-1,     param["r"]: [0.  0.5 1. ]
        geometry = NumbaRadial1DGeometry(
            r[:-1],
            r[1:],
            r[:-1] * c.c.cgs.value,
            r[1:] * c.c.cgs.value,
        )
        return geometry

    @property
    def time_explosion(self):
        # previously used model value that passes tests
        # time taken for a photon to move 1 cm
        return 1 / c.c.cgs.value

    @parameterize({"p": [0.0, 0.5, 1.0], "Test data": TESTDATA})
    def time_calculate_z(self, p, test_data):
        func = formal_integral.calculate_z
        inv_t = 1.0 / self.time_explosion
        len(self.formal_integral_geometry(test_data).r_outer)
        r_outer = self.formal_integral_geometry(test_data).r_outer

        for r in r_outer:
            actual = func(r, p, inv_t)
            if p >= r:
                assert actual == 0
            else:
                np.sqrt(r * r - p * p) * formal_integral.C_INV * inv_t

    @skip_benchmark
    @parameterize({"p": [0, 0.5, 1], "Test data": TESTDATA})
    def time_populate_z_photosphere(self, p, test_data):
        formal_integral.FormalIntegrator(
            self.formal_integral_geometry(test_data), None, None
        )
        func = formal_integral.populate_z
        size = len(self.formal_integral_geometry(test_data).r_outer)
        r_inner = self.formal_integral_geometry(test_data).r_inner
        self.formal_integral_geometry(test_data).r_outer

        p = r_inner[0] * p
        oz = np.zeros_like(r_inner)
        oshell_id = np.zeros_like(oz, dtype=np.int64)

        n = func(
            self.formal_integral_geometry(test_data),
            self.formal_integral_geometry(test_data),
            p,
            oz,
            oshell_id,
        )
        assert n == size

    @skip_benchmark
    @parameterize({"p": [1e-5, 0.5, 0.99, 1], "Test data": TESTDATA})
    def time_populate_z_shells(self, p, test_data):
        formal_integral.FormalIntegrator(
            self.formal_integral_geometry(test_data), None, None
        )
        func = formal_integral.populate_z

        size = len(self.formal_integral_geometry(test_data).r_inner)
        r_inner = self.formal_integral_geometry(test_data).r_inner
        r_outer = self.formal_integral_geometry(test_data).r_outer

        p = r_inner[0] + (r_outer[-1] - r_inner[0]) * p
        idx = np.searchsorted(r_outer, p, side="right")

        oz = np.zeros(size * 2)
        oshell_id = np.zeros_like(oz, dtype=np.int64)

        offset = size - idx

        expected_n = (offset) * 2
        expected_oz = np.zeros_like(oz)
        expected_oshell_id = np.zeros_like(oshell_id)

        # Calculated way to determine which shells get hit
        expected_oshell_id[:expected_n] = (
            np.abs(np.arange(0.5, expected_n, 1) - offset) - 0.5 + idx
        )

        expected_oz[0:offset] = 1 + self.calculate_z(
            r_outer[np.arange(size, idx, -1) - 1], p
        )
        expected_oz[offset:expected_n] = 1 - self.calculate_z(
            r_outer[np.arange(idx, size, 1)], p
        )

        n = func(
            self.formal_integral_geometry(test_data),
            self.formal_integral_geometry(test_data),
            p,
            oz,
            oshell_id,
        )

        assert n == expected_n

    @parameterize({"N": [100, 1000, 10000]})
    def time_calculate_p_values(self, n):
        r = 1.0
        func = formal_integral.calculate_p_values

        expected = r / (n - 1) * np.arange(0, n, dtype=np.float64)
        np.zeros_like(expected, dtype=np.float64)

        func(r, n)
