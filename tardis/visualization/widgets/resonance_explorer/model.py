"""Exact, illustrative inputs for a single-shell resonance flight."""

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class ExplorerInputs:
    """One shell, one directed ray, and one projected-velocity target.

    Lengths and velocities use consistent illustrative units.
    """

    r_inner: float
    r_outer: float
    v_inner: float
    v_outer: float
    impact_parameter: float
    launch_position: float
    line_velocity: float
    s_stop: float | None = None
    stop_cause: str = "assumed electron scatter"

    @property
    def delta_r(self) -> float:
        """Return shell thickness in length units."""
        return self.r_outer - self.r_inner

    @property
    def delta_v(self) -> float:
        """Return boundary velocity difference in velocity units."""
        return self.v_outer - self.v_inner

    @property
    def velocity_slope(self) -> float:
        """Return the radial velocity slope."""
        return self.delta_v / self.delta_r

    @property
    def velocity_intercept(self) -> float:
        """Return the radial intercept using packet transport operation order."""
        # The equivalent inner-edge subtraction can round differently and
        # select another exact Q=0 solver branch for decimal edge values.
        return self.v_outer - self.velocity_slope * self.r_outer

    @property
    def r_0(self) -> float:
        """Return the launch radius."""
        return math.hypot(self.launch_position, self.impact_parameter)

    @property
    def mu_0(self) -> float:
        """Return the launch direction cosine."""
        return self.launch_position / self.r_0
