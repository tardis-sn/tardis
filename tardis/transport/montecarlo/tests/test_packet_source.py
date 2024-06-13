import os

from astropy import units as u
import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from tardis.transport.montecarlo.packet_source import (
    BlackBodySimpleSource,
    BlackBodySimpleSourceRelativistic,
)
from tardis.transport.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
from tardis.tests.fixtures.regression_data import RegressionData


class TestBlackBodySimpleSource:
    @pytest.fixture(scope="class")
    def blackbodysimplesource(self, request):
        """
        Create BlackBodySimpleSource instance.

        Yields
        -------
        tardis.transport.montecarlo.packet_source.BlackBodySimpleSource
        """
        cls = type(self)
        bb = BlackBodySimpleSource(
            radius=123,
            temperature=10000 * u.K,
            base_seed=1963,
            legacy_second_seed=2508,
            legacy_mode_enabled=True,
        )
        yield bb

    def test_bb_nus(self, regression_data, blackbodysimplesource):
        actual_nus = blackbodysimplesource.create_packet_nus(100).value
        expected_nus = regression_data.sync_ndarray(actual_nus)
        assert_allclose(actual_nus, expected_nus)

    def test_bb_mus(self, regression_data, blackbodysimplesource):
        actual_mus = blackbodysimplesource.create_packet_mus(100)
        expected_mus = regression_data.sync_ndarray(actual_mus)
        assert_allclose(actual_mus, expected_mus)

    def test_bb_energies(self, regression_data, blackbodysimplesource):
        actual_unif_energies = blackbodysimplesource.create_packet_energies(
            100
        ).value
        expected_unif_energies = regression_data.sync_ndarray(
            actual_unif_energies
        )
        assert_allclose(actual_unif_energies, expected_unif_energies)

    def test_bb_attributes(self, regression_data, blackbodysimplesource):
        actual_bb = blackbodysimplesource
        expected_bb = regression_data.sync_hdf_store(actual_bb)[
            "/black_body_simple_source/scalars"
        ]
        assert_allclose(expected_bb.base_seed, actual_bb.base_seed)
        assert_allclose(expected_bb.temperature, actual_bb.temperature.value)
        assert_allclose(expected_bb.radius, actual_bb.radius)


class TestBlackBodySimpleSourceRel:
    @pytest.fixture(scope="class")
    def blackbody_simplesource_relativistic(self):
        """
        Create BlackBodySimpleSourceRelativistic instance.

        Yields
        -------
        tardis.montecarlo.packet_source.BlackBodySimpleSourceRelativistic
        """
        bb_rel = BlackBodySimpleSourceRelativistic(
            time_explosion=1123187,
            base_seed=1963,
            legacy_second_seed=2508,
            legacy_mode_enabled=True,
        )
        bb_rel.temperature = 10000 * u.K
        bb_rel.beta = 0.25
        yield bb_rel

    def test_bb_nus(self, regression_data, blackbody_simplesource_relativistic):
        actual_nus = blackbody_simplesource_relativistic.create_packet_nus(
            100
        ).value
        expected_nus = regression_data.sync_ndarray(actual_nus)
        assert_allclose(actual_nus, expected_nus)

    def test_bb_energies(
        self, regression_data, blackbody_simplesource_relativistic
    ):
        actual_unif_energies = (
            blackbody_simplesource_relativistic.create_packet_energies(
                100
            ).value
        )
        gamma = np.sqrt(1 - blackbody_simplesource_relativistic.beta**2) ** -1
        expected_unif_energies = regression_data.sync_ndarray(
            actual_unif_energies
        )
        assert_allclose(actual_unif_energies, expected_unif_energies)

    def test_bb_mus(self, regression_data, blackbody_simplesource_relativistic):
        blackbody_simplesource_relativistic._reseed(2508)
        actual_mus = blackbody_simplesource_relativistic.create_packet_mus(10)
        expected_mus = regression_data.sync_ndarray(actual_mus)
        assert_allclose(actual_mus, expected_mus)

    def test_bb_attributes(
        self, regression_data, blackbody_simplesource_relativistic
    ):
        actual_bb = blackbody_simplesource_relativistic
        expected_bb = regression_data.sync_hdf_store(actual_bb)[
            "/black_body_simple_source/scalars"
        ]
        assert_allclose(expected_bb.base_seed, actual_bb.base_seed)
        assert_allclose(expected_bb.temperature, actual_bb.temperature.value)
        assert_allclose(expected_bb.time_explosion, actual_bb.time_explosion)
