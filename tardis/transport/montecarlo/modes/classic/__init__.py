"""Classic Monte Carlo transport mode.

Classic mode implements line-only radiation transport without continuum processes.
This is the standard TARDIS transport mode used for most supernova simulations.
"""

from tardis.transport.montecarlo.modes.classic.montecarlo_transport import (
    montecarlo_transport,
)

__all__ = ["montecarlo_transport"]
