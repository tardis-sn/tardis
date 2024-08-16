from tardis.opacities.tau_sobolev import calculate_sobolev_line_opacity
from tardis.opacities.opacity_state import (
    OpacityState,
)
import numpy as np
import pandas as pd


class OpacitySolver(object):

    line_interaction_type: str = "scatter"
    disable_line_scattering: bool = False

    def __init__(
        self, line_interaction_type="scatter", disable_line_scattering=False
    ):
        """Solver class for opacities

        Parameters
        ----------
        line_interaction_type: str
            "scatter", "downbranch", or "macroatom"
        disable_line_scattering: bool
        """

        self.line_interaction_type = line_interaction_type
        self.disable_line_scattering = disable_line_scattering

    def solve(self, legacy_plasma) -> OpacityState:
        """
        Solves the opacity state

        Parameters
        ----------
        plasma : tarids.plasma.BasePlasma
            legacy base plasma

        Returns
        -------
        OpacityState
        """
        atomic_data = legacy_plasma.atomic_data

        if self.disable_line_scattering:
            tau_sobolev = pd.DataFrame(
                np.zeros(
                    (
                        legacy_plasma.atomic_data.lines.shape[
                            0
                        ],  # number of lines
                        legacy_plasma.number_density.shape[
                            1
                        ],  # number of shells
                    ),
                    dtype=np.float64,
                ),
                index=atomic_data.lines.index,
            )
        else:
            tau_sobolev = calculate_sobolev_line_opacity(
                atomic_data.lines,
                legacy_plasma.level_number_density,
                legacy_plasma.time_explosion,
                legacy_plasma.stimulated_emission_factor,
            )

        opacity_state = OpacityState.from_legacy_plasma(
            legacy_plasma, tau_sobolev
        )

        return opacity_state
