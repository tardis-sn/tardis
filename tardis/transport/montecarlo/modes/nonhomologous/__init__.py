"""Non-homologous Monte Carlo transport mode.

Non-homologous mode is based on classic mode, which implements line-only radiation
transport without continuum processes. Non-homologous mode includes more generalized
geometry and transport to handle the case where the velocity profile of the ejecta
is a non-monotonic function of radius.
"""

from tardis.transport.montecarlo.modes.nonhomologous.montecarlo_transport import (
    montecarlo_transport,
)

__all__ = ["montecarlo_transport"]
