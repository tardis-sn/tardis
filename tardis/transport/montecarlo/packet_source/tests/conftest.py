"""Configuration for packet source tests."""

import pytest
from astropy import units as u

from tardis.transport.montecarlo.packet_source.black_body_weighted import (
    BlackBodyWeightedSource,
)

# Import all fixtures from the parent montecarlo tests conftest
# This ensures the packet source tests have access to all necessary fixtures
from tardis.transport.montecarlo.tests.conftest import *  # noqa: F403


@pytest.fixture(scope="session")
def simple_weighted_packet_source():
    """Create a simple weighted packet source for testing."""
    return BlackBodyWeightedSource(
        radius=1e14 * u.cm,
        temperature=10000 * u.K,
        base_seed=1963,
        legacy_mode_enabled=True,
    )
