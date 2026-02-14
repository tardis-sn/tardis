"""Test IIP Monte Carlo transport mode."""

import pytest


def test_iip_mode_imports():
    """Test that IIP mode can be imported correctly."""
    from tardis.transport.montecarlo.modes.iip import montecarlo_transport

    assert callable(montecarlo_transport), (
        "montecarlo_transport should be callable"
    )


def test_iip_mode_returns_5_tuple(simulation_verysimple):
    """
    Test that IIP mode returns 5-tuple including continuum estimators.

    This test requires continuum to be enabled in the simulation config.
    """
    # This is a placeholder - actual implementation would:
    # 1. Enable continuum in simulation
    # 2. Run IIP mode transport
    # 3. Verify 5-tuple return
    # 4. Verify estimators_continuum is populated
    pytest.skip("Requires continuum-enabled simulation fixture")
