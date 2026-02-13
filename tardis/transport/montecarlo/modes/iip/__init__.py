"""IIP Monte Carlo transport mode.

IIP mode implements full radiation transport with both continuum processes and
full relativistic corrections ALWAYS ENABLED. This mode is optimized for Type IIP
supernova simulations requiring accurate continuum interactions and relativistic effects.

Key Features:
- Continuum processes (bound-free, free-free, free-bound) always active
- Full relativistic transformations (Doppler shifts, aberration) always enabled
- Returns 5-tuple including continuum estimators
- No conditional logic for continuum or relativity flags

This mode is a specialized, performance-optimized implementation where advanced
physics is mandatory rather than optional.
"""

from .montecarlo_transport import montecarlo_transport

__all__ = ["montecarlo_transport"]
