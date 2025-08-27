"""Configuration for packet source tests."""

import pytest
from astropy import units as u

from tardis.transport.montecarlo.packet_source.weighted import (
    BlackBodyWeightedSource,
)

# Import fixtures from the parent montecarlo tests directory
pytest_plugins = ["tardis.transport.montecarlo.tests.conftest"]


@pytest.fixture(scope="session")
def simple_weighted_packet_source():
    """Create a simple weighted packet source for testing."""
    return BlackBodyWeightedSource(
        radius=1e14 * u.cm,
        temperature=10000 * u.K,
        base_seed=1963,
        legacy_mode_enabled=True,
    )
