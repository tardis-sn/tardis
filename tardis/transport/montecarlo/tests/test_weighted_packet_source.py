
from astropy import units as u
import pytest
from numpy.testing import assert_allclose

from tardis.transport.montecarlo.weighted_packet_source import (
    BlackBodyWeightedSource,
)


class TestBlackBodyWeightedSource:
    @pytest.fixture(scope="class")
    def blackbodyweightedsource(self, request):
        """
        Create blackbodyweightedsource instance.

        Yields
        -------
        tardis.transport.montecarlo.packet_source.blackbodyweightedsource
        """
        cls = type(self)
        bb = BlackBodyWeightedSource(
            radius=123,
            temperature=10000 * u.K,
            base_seed=1963,
            legacy_mode_enabled=False,
        )
        yield bb

    def test_bb_nus(self, regression_data, blackbodyweightedsource):
        actual_nus = blackbodyweightedsource.create_packet_nus(100).value
        expected_nus = regression_data.sync_ndarray(actual_nus)
        assert_allclose(actual_nus, expected_nus)

    def test_bb_mus(self, regression_data, blackbodyweightedsource):
        actual_mus = blackbodyweightedsource.create_packet_mus(100)
        expected_mus = regression_data.sync_ndarray(actual_mus)
        assert_allclose(actual_mus, expected_mus)

    def test_bb_energies(self, regression_data, blackbodyweightedsource):
        actual_unif_energies = blackbodyweightedsource.create_packet_energies(
            100
        ).value
        expected_unif_energies = regression_data.sync_ndarray(
            actual_unif_energies
        )
        assert_allclose(actual_unif_energies, expected_unif_energies)

    def test_bb_attributes(self, regression_data, blackbodyweightedsource):
        actual_bb = blackbodyweightedsource
        expected_bb = regression_data.sync_hdf_store(actual_bb)[
            "/black_body_weighted_source/scalars"
        ]
        assert_allclose(expected_bb.base_seed, actual_bb.base_seed)
        assert_allclose(expected_bb.temperature, actual_bb.temperature.value)
        assert_allclose(expected_bb.radius, actual_bb.radius)
