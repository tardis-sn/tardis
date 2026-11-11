"""Exact, illustrative inputs for a single-shell resonance flight."""

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class ExplorerInputs:
    """One shell, one directed ray, and one projected-velocity target.

    Lengths and velocities use consistent illustrative units.
    """

    r_in: float
    r_out: float
    v_in: float
    v_out: float
    b: float
    x_0: float
    w: float
    s_stop: float | None = None
    stop_cause: str = "assumed electron scatter"

    @property
    def delta_r(self) -> float:
        """Return shell thickness in length units."""
        return self.r_out - self.r_in

    @property
    def delta_v(self) -> float:
        """Return boundary velocity difference in velocity units."""
        return self.v_out - self.v_in

    @property
    def m(self) -> float:
        """Return the radial velocity slope."""
        return self.delta_v / self.delta_r

    @property
    def q(self) -> float:
        """Return the radial intercept using packet transport operation order."""
        # The equivalent inner-edge subtraction can round differently and
        # select another exact Q=0 solver branch for decimal edge values.
        return self.v_out - self.m * self.r_out

    @property
    def r_0(self) -> float:
        """Return the launch radius."""
        return math.hypot(self.x_0, self.b)

    @property
    def mu_0(self) -> float:
        """Return the launch direction cosine."""
        return self.x_0 / self.r_0
