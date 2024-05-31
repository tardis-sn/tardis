import numpy as np
import pytest

from tardis.energy_input.gamma_ray_packet_source import RadioactivePacketSource


@pytest.mark.skip(reason="Packet source init is very complex")
class TestGammaRayPacketSource:
    @pytest.fixture(scope="class")
    def radioactivepacketsource(self, request):
        """
        Create RadioactivePacketSource instance.

        Yields
        -------
        tardis.energy_input.gamma_ray_packet_source.RadioactivePacketSource
        """
        cls = type(self)
        cls.packet_source = RadioactivePacketSource(base_seed=1963)
        yield cls.packet_source

    def test_create_packet_radii(
        self, regression_data, radioactivepacketsource
    ):
        actual = self.packet_source.create_packet_radii()
        expected = regression_data.sync_ndarray(actual)
        assert np.all(np.isclose(actual, expected))

    def test_create_packet_nus(self, regression_data, radioactivepacketsource):
        actual = self.packet_source.create_packet_nus()
        expected = regression_data.sync_ndarray(actual)
        assert np.all(np.isclose(actual, expected))

    def test_create_packet_directions(
        self, regression_data, radioactivepacketsource
    ):
        actual = self.packet_source.create_packet_directions()
        expected = regression_data.sync_ndarray(actual)
        assert np.all(np.isclose(actual, expected))

    def test_create_packet_energies(
        self, regression_data, radioactivepacketsource
    ):
        actual = self.packet_source.create_packet_energies()
        expected = regression_data.sync_ndarray(actual)
        assert np.all(np.isclose(actual, expected))

    def test_create_packet_times_uniform_time(
        self, regression_data, radioactivepacketsource
    ):
        actual = self.packet_source.create_packet_times_uniform_time()
        expected = regression_data.sync_ndarray(actual)
        assert np.all(np.isclose(actual, expected))

    def test_create_packet_times_uniform_energy(
        self, regression_data, radioactivepacketsource
    ):
        actual = self.packet_source.create_packet_times_uniform_energy()
        expected = regression_data.sync_ndarray(actual)
        assert np.all(np.isclose(actual, expected))

    def test_calculate_energy_factors(
        self, regression_data, radioactivepacketsource
    ):
        actual = self.packet_source.calculate_energy_factors()
        expected = regression_data.sync_ndarray(actual)
        assert np.all(np.isclose(actual, expected))

    def test_create_packets(self, regression_data, radioactivepacketsource):
        assert True

    def test_calculate_positron_fraction(
        self, regression_data, radioactivepacketsource
    ):
        actual = self.packet_source.calculate_positron_fraction()
        expected = regression_data.sync_ndarray(actual)
        assert np.all(np.isclose(actual, expected))
