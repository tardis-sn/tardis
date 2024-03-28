"""
Basic TARDIS Benchmark.
"""
from asv_runner.benchmarks.mark import skip_benchmark, parameterize
from numpy.testing import (
    assert_almost_equal,
)

import tardis.transport.frame_transformations as frame_transformations
from benchmarks.benchmark_base import BenchmarkBase


# @skip_benchmark
class BenchmarkTransportDopplerFactor(BenchmarkBase):
    """
    Class to benchmark the Doppler factor function.
    """

    def __init__(self):
        pass

    @staticmethod
    @parameterize({"Parameters": [
        {
            "mu": 0.3,
            "r": 7.5e14,
            "inv_t_exp": 1 / 5.2e7,
            "expected": 0.9998556693818854
        },
        {
            "mu": -0.3,
            "r": 0,
            "inv_t_exp": 1 / 2.6e7,
            "expected": 1.0
        },
        {
            "mu": 0,
            "r": 1,
            "inv_t_exp": 1 / 2.6e7,
            "expected": 1.0
        },
    ]})
    def time_get_doppler_factor(parameters):
        mu = parameters["mu"]
        r = parameters["r"]
        inv_t_exp = parameters["inv_t_exp"]
        expected = parameters["expected"]
        # Set the params from test cases here
        time_explosion = 1 / inv_t_exp

        # Perform any other setups just before this, they can be additional calls
        # to other methods or introduction of some temporary variables

        obtained = frame_transformations.get_doppler_factor(r, mu, time_explosion)

        # Perform required assertions
        assert_almost_equal(obtained, expected)

    @staticmethod
    @parameterize({"Parameters": [
        {
            "mu": 0.3,
            "beta": 0.2,
            "expected": 0.94
        },
        {
            "mu": -0.3,
            "beta": 0,
            "expected": 1.0
        },
        {
            "mu": 0,
            "beta": 0.8,
            "expected": 1.0
        },
    ]})
    def time_get_doppler_factor_partial_relativity(parameters):
        mu = parameters["mu"]
        beta = parameters["beta"]
        expected = parameters["expected"]
        obtained = frame_transformations.get_doppler_factor_partial_relativity(
            mu, beta
        )
        assert_almost_equal(obtained, expected)

    @staticmethod
    @parameterize({"Parameters": [
        {
            "mu": 0.3,
            "beta": 0.2,
            "expected": 0.95938348
        },
        {
            "mu": -0.3,
            "beta": 0,
            "expected": 1.0
        },
        {
            "mu": 0,
            "beta": 0.8,
            "expected": 1.6666667
        },
    ]})
    def time_get_doppler_factor_full_relativity(parameters):
        mu = parameters["mu"]
        beta = parameters["beta"]
        expected = parameters["expected"]
        obtained = frame_transformations.get_doppler_factor_full_relativity(
            mu, beta
        )
        assert_almost_equal(obtained, expected)

    @staticmethod
    @parameterize({"Parameters": [
        {
            "mu": 0.3,
            "r": 7.5e14,
            "inv_t_exp": 1 / 5.2e7,
            "expected": 1 / 0.9998556693818854
        },
        {
            "mu": -0.3,
            "r": 0,
            "inv_t_exp": 1 / 2.6e7,
            "expected": 1.0
        },
        {
            "mu": 0,
            "r": 1,
            "inv_t_exp": 1 / 2.6e7,
            "expected": 1.0
        },
    ]})
    def time_get_inverse_doppler_factor(parameters):
        mu = parameters["mu"]
        r = parameters["r"]
        inv_t_exp = parameters["inv_t_exp"]
        expected = parameters["expected"]
        # Set the params from test cases here
        time_explosion = 1 / inv_t_exp

        # Perform any other setups just before this, they can be additional calls
        # to other methods or introduction of some temporary variables

        obtained = frame_transformations.get_inverse_doppler_factor(
            r, mu, time_explosion
        )

        # Perform required assertions
        assert_almost_equal(obtained, expected)

    @staticmethod
    @skip_benchmark
    @parameterize({"Parameters": [
        ["mu", "beta", "expected"],
        [
            (0.3, 0.2, 1 / 0.94),
            (-0.3, 0, 1.0),
            (0, 0.8, 1.0),
        ],
    ]})
    def time_get_inverse_doppler_factor_partial_relativity(parameters):
        mu = parameters["mu"]
        beta = parameters["beta"]
        expected = parameters["expected"]
        obtained = (
            frame_transformations.get_inverse_doppler_factor_partial_relativity(
                mu, beta
            )
        )
        assert_almost_equal(obtained, expected)

    @staticmethod
    @parameterize({"Parameters": [
        {
            "mu": 0.3,
            "beta": 0.2,
            "expected": 1.0818579
        },
        {
            "mu": -0.3,
            "beta": 0,
            "expected": 1.0
        },
        {
            "mu": 0,
            "beta": 0.8,
            "expected": 1.6666667
        },
    ]})
    def time_get_inverse_doppler_factor_full_relativity(parameters):
        mu = parameters["mu"]
        beta = parameters["beta"]
        expected = parameters["expected"]
        obtained = frame_transformations.get_inverse_doppler_factor_full_relativity(
            mu, beta
        )
        assert_almost_equal(obtained, expected)
