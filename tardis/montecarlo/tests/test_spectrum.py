import os
import pytest
import numpy as np
import pandas as pd

from astropy import (
        units as u,
        constants as c
        )
import astropy.tests.helper as test_helper
from tardis.montecarlo.spectrum import TARDISSpectrum

BIN = 5

TESTDATA = [
        {
            'nu': u.Quantity(np.linspace(1, 50, BIN + 1), 'Hz'),
            'lum':  u.Quantity(np.linspace(1e27, 1e28, BIN), 'erg / s'),
            'distance': None
            },
        {
            'nu': u.Quantity(np.linspace(1, 50, BIN + 1), 'Hz'),
            'lum':  u.Quantity(np.linspace(1e27, 1e28, BIN), 'erg / s'),
            'distance': u.Quantity(1, 'Mpc')
            },
        {
            'nu': u.Quantity([1, 2, 3, 4], 'Hz'),
            'lum': u.Quantity([1, 2, 3], 'erg / s') * np.pi,
            'distance': u.Quantity(.5, 'cm'),
            }
        ]


@pytest.fixture(
        scope='module',
        params=TESTDATA
        )
def spectrum(request):
    data = TARDISSpectrum(
            request.param['nu'],
            request.param['lum'],
            )
    distance = request.param['distance']
    if distance:
        data.distance = distance
    return data


###
# Testing properties
###
def test_frequency(spectrum):
    assert np.all(
            spectrum.frequency == spectrum._frequency[:-1]
            )


def test_luminosity_density_nu(spectrum):
    expected = spectrum.luminosity / np.diff(spectrum._frequency)
    test_helper.assert_quantity_allclose(
            spectrum.luminosity_density_nu,
            expected
            )


def test_luminosity_density_lambda(spectrum):
    expected = spectrum.f_nu_to_f_lambda(
            spectrum.luminosity / np.diff(spectrum._frequency)
            )
    test_helper.assert_quantity_allclose(
            spectrum.luminosity_density_lambda,
            expected
            )


def test_flux_nu(spectrum):
    if getattr(spectrum, 'distance', None):
        test_helper.assert_quantity_allclose(
                spectrum.flux_nu,
                spectrum.luminosity_to_flux(
                    spectrum.luminosity_density_nu,
                    spectrum.distance)
                )
    else:
        with pytest.raises(AttributeError):
            spectrum.flux_nu


def test_flux_lambda(spectrum):
    if getattr(spectrum, 'distance', None):
        test_helper.assert_quantity_allclose(
                spectrum.flux_lambda,
                spectrum.luminosity_to_flux(
                    spectrum.luminosity_density_lambda,
                    spectrum.distance)
                )
    else:
        with pytest.raises(AttributeError):
            spectrum.flux_nu


###
# Testing methods
###

def test_luminosity_to_flux():
    lum = u.Quantity(np.arange(1, 4, 1) * np.pi, 'erg / s')
    distance = u.Quantity(.5, 'cm')
    flux = TARDISSpectrum.luminosity_to_flux(lum, distance)
    test_helper.assert_quantity_allclose(
            flux,
            u.Quantity(lum.value / np.pi, 'erg/s/cm^2')
            )


def test_f_nu_to_f_lambda(spectrum):
    expected = (
            spectrum.luminosity_density_nu *
            spectrum.frequency**2 / c.c.to('angstrom/s')
            ).to('erg / (s angstrom)')
    print(expected)
    test_helper.assert_quantity_allclose(
            spectrum.luminosity_density_lambda,
            expected
            )
    expected = (
            spectrum.luminosity_density_nu *
            spectrum.frequency**2 / c.c.to('angstrom/s')
            ).to('erg / (s angstrom)')
    np.testing.assert_allclose(
            spectrum.luminosity_density_lambda.value,
            expected.value
            )


###
# Save and Load
###

def compare_spectra(actual, desired):
    test_helper.assert_quantity_allclose(
            actual.frequency,
            desired.frequency
            )
    test_helper.assert_quantity_allclose(
            actual.luminosity,
            desired.luminosity
            )
    if getattr(actual, 'distance', None):
        test_helper.assert_quantity_allclose(
                actual.distance,
                desired.distance
                )


@pytest.mark.xfail
def test_to_from_hdf(tmpdir, spectrum):
    path = str(tmpdir.join('spectrum.hdf'))
    spectrum.to_hdf(
            str(path),
            'spectrum'
            )

    spec_from = TARDISSpectrum.from_hdf(
            path,
            'spectrum'
            )

    compare_spectra(
            spectrum,
            spec_from)


@pytest.mark.xfail
def test_to_from_hdf_buffer(tmpdir, spectrum):
    path = str(tmpdir.join('spectrum.hdf'))
    spectrum.to_hdf(
            str(path),
            'spectrum'
            )

    with pd.HDFStore(path, mode='r') as buffer:
        spec_from = TARDISSpectrum.from_hdf(
                buffer,
                'spectrum'
                )

    compare_spectra(
            spectrum,
            spec_from)

def test_creat_from_wl(spectrum):
    actual = TARDISSpectrum(
            spectrum._frequency.to('angstrom', u.spectral()),
            spectrum.luminosity
            )

    compare_spectra(actual, spectrum)
