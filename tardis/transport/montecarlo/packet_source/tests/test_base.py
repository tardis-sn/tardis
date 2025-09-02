import numpy as np
import pytest
from astropy import units as u

from tardis.transport.montecarlo.packet_source.base import BasePacketSource


class TestablePacketSource(BasePacketSource):
    """Concrete implementation of BasePacketSource for testing."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.radius = 1e14 * u.cm
        self.temperature = 10000 * u.K

    def create_packet_nus(self, no_of_packets, *args, **kwargs):
        return np.ones(no_of_packets) * 1e15 * u.Hz

    def create_packet_mus(self, no_of_packets, *args, **kwargs):
        return np.ones(no_of_packets) * 0.5

    def create_packet_energies(self, no_of_packets, *args, **kwargs):
        return np.ones(no_of_packets) * u.erg

    def create_packet_radii(self, no_of_packets, *args, **kwargs):
        return np.ones(no_of_packets) * self.radius


class TestBasePacketSource:
    @pytest.fixture
    def base_packet_source(self):
        """Create a testable packet source instance."""
        return TestablePacketSource(base_seed=1963)

    def test_initialization(self, base_packet_source):
        """Test basic initialization of packet source."""
        assert base_packet_source.base_seed == 1963
        assert base_packet_source.MAX_SEED_VAL == 2**32 - 1

    def test_reseed(self, base_packet_source):
        """Test reseeding functionality."""
        base_packet_source._reseed(42)
        # Check that rng is created and is different from initial state
        assert hasattr(base_packet_source, 'rng')
        assert base_packet_source.rng is not None

    def test_create_packets(self, base_packet_source):
        """Test packet creation."""
        packets = base_packet_source.create_packets(10)

        # Check that packet collection is created
        assert packets is not None
        assert len(packets.initial_radii) == 10
        assert len(packets.initial_nus) == 10
        assert len(packets.initial_mus) == 10
        assert len(packets.initial_energies) == 10
        assert len(packets.packet_seeds) == 10

        # Check that all values are reasonable
        assert np.all(packets.initial_radii > 0)
        assert np.all(packets.initial_nus > 0)
        assert np.all(packets.initial_energies > 0)
        assert np.all(np.abs(packets.initial_mus) <= 1)

    def test_create_packets_no_base_seed(self):
        """Test that packet creation fails without base_seed."""
        source = TestablePacketSource(base_seed=None)
        with pytest.raises(ValueError, match="base_seed must be set"):
            source.create_packets(10)

    def test_calculate_radfield_luminosity(self, base_packet_source):
        """Test radiation field luminosity calculation."""
        luminosity = base_packet_source.calculate_radfield_luminosity()

        # Check that luminosity is positive and has correct units
        assert luminosity > 0
        assert luminosity.unit.is_equivalent(u.erg / u.s)

        # Check that luminosity scales with temperature^4 and radius^2
        base_packet_source.temperature *= 2
        luminosity_2x_temp = base_packet_source.calculate_radfield_luminosity()
        assert (luminosity_2x_temp / luminosity).value == pytest.approx(16, rel=1e-10)

    def test_abstract_methods_must_be_implemented(self):
        """Test that abstract methods must be implemented in subclasses."""
        with pytest.raises(TypeError):
            BasePacketSource()  # Cannot instantiate abstract class
