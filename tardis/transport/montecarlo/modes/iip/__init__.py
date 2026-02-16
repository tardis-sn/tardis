"""
IIP Monte Carlo Transport Mode.

This module implements the Monte Carlo transport for Type IIP supernovae,
with continuum processes and full relativistic corrections always enabled.
"""

from tardis.transport.montecarlo.modes.iip.montecarlo_transport import (
    montecarlo_transport as montecarlo_transport,
)

__all__ = ["montecarlo_transport"]
