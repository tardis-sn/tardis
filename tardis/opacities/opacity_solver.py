from tardis.opacities.tau_sobolev import calculate_sobolev_line_opacity
from tardis.opacities.opacity_state import (
    OpacityStatePython,
)
import numpy as np


class OpacitySolver(object):
    def __init__(self, line_interaction_type="scatter", disable_line_scattering=False):
        """Solver class for opacities

        Parameters
        ----------
        line_interaction_type: str
            "scatter", "downbranch", or "macroatom"
        disable_line_scattering: bool
        """

        self.line_interaction_type = line_interaction_type
        self.disable_line_scattering = disable_line_scattering

    def solve(self, legacy_plasma) -> OpacityStatePython:
        """
        Solves the opacity state

        Parameters
        ----------
        plasma : tarids.plasma.BasePlasma
            legacy base plasma

        Returns
        -------
        OpacityStatePython
        """
        tau_sobolev = calculate_sobolev_line_opacity(
            legacy_plasma.atomic_data.lines,
            legacy_plasma.level_number_density,
            legacy_plasma.time_explosion,
            legacy_plasma.stimulated_emission_factor,
        )

        if self.disable_line_scattering:
            tau_sobolev *= 0.0

        opacity_state = OpacityStatePython.from_legacy_plasma(
            legacy_plasma, tau_sobolev
        )

        return opacity_state
