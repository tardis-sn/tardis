"""
Basic TARDIS Benchmark.
"""
import os

import astropy.tests.helper as test_helper
import numpy as np
import pandas as pd
from astropy import units as u
from asv_runner.benchmarks.mark import parameterize, skip_benchmark
from numpy.testing import assert_almost_equal

from benchmarks.benchmark_base import BenchmarkBase
from tardis import constants as c
from tardis.montecarlo.spectrum import (
    TARDISSpectrum,
)


# @skip_benchmark
class BenchmarkMontecarloSpectrum(BenchmarkBase):
    """
    Class to benchmark the spectrum function.
    """

    def __init__(self):
        pass

    BIN = 5

    TESTDATA = [
        {
            "nu": u.Quantity(np.linspace(1, 50, BIN + 1), "Hz"),
            "lum": u.Quantity(np.linspace(1e27, 1e28, BIN), "erg / s"),
            "distance": None,
        },
        {
            "nu": u.Quantity(np.linspace(1, 50, BIN + 1), "Hz"),
            "lum": u.Quantity(np.linspace(1e27, 1e28, BIN), "erg / s"),
            "distance": u.Quantity(1, "Mpc"),
        },
        {
            "nu": u.Quantity([1, 2, 3, 4], "Hz"),
            "lum": u.Quantity([1, 2, 3], "erg / s") * np.pi,
            "distance": u.Quantity(0.5, "cm"),
        },
    ]

    def spectrum(self, param):
        data = TARDISSpectrum(
            param["nu"],
            param["lum"],
        )
        distance = param["distance"]
        if distance is not None:
            data.distance = distance
        return data

    ###
    # Testing properties
    ###

    @parameterize({"Spectrum": TESTDATA})
    def time_frequency(self, parameters):
        spectrum = self.spectrum(parameters)
        assert np.all(spectrum.frequency == spectrum._frequency[:-1])

    @parameterize({"Spectrum": TESTDATA})
    def time_luminosity_density_nu(self, parameters):
        spectrum = self.spectrum(parameters)
        expected = spectrum.luminosity / np.diff(spectrum._frequency)
        test_helper.assert_quantity_allclose(
            spectrum.luminosity_density_nu, expected
        )

    @parameterize({"Spectrum": TESTDATA})
    def time_luminosity_density_lambda(self, parameters):
        spectrum = self.spectrum(parameters)
        expected = spectrum.f_nu_to_f_lambda(
            spectrum.luminosity / np.diff(spectrum._frequency)
        )
        test_helper.assert_quantity_allclose(
            spectrum.luminosity_density_lambda, expected
        )

    # TODO: Make decisions about how to implement the PyTest Warns in benchmarks.
    # @parameterize({"Spectrum": TESTDATA})
    # def time_flux_nu(self, parameters):
    #     spectrum = self.spectrum(parameters)
    #     if getattr(spectrum, "distance", None) is not None:
    #
    #         with pytest.warns(DeprecationWarning):
    #             test_helper.assert_quantity_allclose(
    #                 spectrum.flux_nu,
    #                 spectrum.luminosity_to_flux(
    #                     spectrum.luminosity_density_nu, spectrum.distance
    #                 ),
    #             )
    #     else:
    #         with pytest.raises(AttributeError):
    #             spectrum.flux_nu

    # TODO: Make decisions about how to implement the PyTest Warns in benchmarks.
    # @parameterize({"Spectrum": TESTDATA})
    # def time_flux_lambda(self, parameters):
    #     spectrum = self.spectrum(parameters)
    #     if getattr(spectrum, "distance", None) is not None:
    #
    #         with pytest.warns(DeprecationWarning):
    #             test_helper.assert_quantity_allclose(
    #                 spectrum.flux_lambda,
    #                 spectrum.luminosity_to_flux(
    #                     spectrum.luminosity_density_lambda, spectrum.distance
    #                 ),
    #             )
    #     else:
    #         with pytest.raises(AttributeError):
    #             spectrum.flux_nu

    ###
    # Testing methods
    ###

    def time_luminosity_to_flux(self):
        lum = u.Quantity(np.arange(1, 4, 1) * np.pi, "erg / s")
        distance = u.Quantity(0.5, "cm")
        flux = TARDISSpectrum.luminosity_to_flux(lum, distance)
        test_helper.assert_quantity_allclose(
            flux, u.Quantity(lum.value / np.pi, "erg s^-1 cm^-2")
        )

    @parameterize({"Spectrum": TESTDATA})
    def time_f_nu_to_f_lambda(self, parameters):
        spectrum = self.spectrum(parameters)
        expected = (
                spectrum.luminosity_density_nu
                * spectrum.frequency ** 2
                / c.c.to("angstrom/s")
        ).to("erg / (s angstrom)")
        print(expected)
        test_helper.assert_quantity_allclose(
            spectrum.luminosity_density_lambda, expected
        )
        expected = (
                spectrum.luminosity_density_nu
                * spectrum.frequency ** 2
                / c.c.to("angstrom/s")
        ).to("erg / (s angstrom)")
        np.testing.assert_allclose(
            spectrum.luminosity_density_lambda.value, expected.value
        )

    ###
    # Save and Load
    ###

    def compare_spectra(self, actual, desired):
        test_helper.assert_quantity_allclose(actual.frequency, desired.frequency)
        test_helper.assert_quantity_allclose(actual.luminosity, desired.luminosity)
        if getattr(actual, "distance", None):
            test_helper.assert_quantity_allclose(actual.distance, desired.distance)

    @parameterize({"Spectrum": TESTDATA, "Attributes": TARDISSpectrum.hdf_properties})
    def time_hdf_spectrum(self, parameters, attr):
        hdf_file_path = self.hdf_file_path
        spectrum = self.spectrum(parameters)
        spectrum.to_hdf(hdf_file_path, name="spectrum", overwrite=True)
        actual = getattr(spectrum, attr)
        if hasattr(actual, "cgs"):
            actual = actual.cgs.value

        if np.isscalar(actual):
            path = os.path.join("spectrum", "scalars")
            expected = getattr(pd.read_hdf(hdf_file_path, path), attr)
            assert_almost_equal(actual, expected)
        else:
            path = os.path.join("spectrum", attr)
            expected = pd.read_hdf(hdf_file_path, path)
            assert_almost_equal(actual, expected.values)

    ###
    # Test creation from nonstandard units
    ###

    @parameterize({"Spectrum": TESTDATA})
    def time_creat_from_wl(self, parameters):
        spectrum = self.spectrum(parameters)
        actual = TARDISSpectrum(
            spectrum._frequency.to("angstrom", u.spectral()), spectrum.luminosity
        )

        self.compare_spectra(actual, spectrum)

    @parameterize({"Spectrum": TESTDATA})
    def time_creat_from_J(self, parameters):
        spectrum = self.spectrum(parameters)
        actual = TARDISSpectrum(
            spectrum._frequency, spectrum.luminosity.to("J / s")
        )

        self.compare_spectra(actual, spectrum)
