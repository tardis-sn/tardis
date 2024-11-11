import numpy as np
import pandas as pd

from tardis.opacities.opacity_state import (
    OpacityState,
)
from tardis.opacities.tau_sobolev import (
    calculate_beta_sobolev,
    calculate_sobolev_line_opacity,
)


class OpacitySolver:
    line_interaction_type: str = "scatter"
    disable_line_scattering: bool = False

    def __init__(
        self,
        line_interaction_type="scatter",
        disable_line_scattering=False,
        previous_beta_sobolev=None,
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
        self.previous_beta_sobolev = previous_beta_sobolev

    def legacy_solve(self, plasma) -> OpacityState:
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
        atomic_data = plasma.atomic_data

        if self.disable_line_scattering:
            tau_sobolev = pd.DataFrame(
                np.zeros(
                    (
                        atomic_data.lines.shape[0],  # number of lines
                        plasma.number_density.shape[1],  # number of shells
                    ),
                    dtype=np.float64,
                ),
                index=atomic_data.lines.index,
            )
        else:
            tau_sobolev = calculate_sobolev_line_opacity(
                atomic_data.lines,
                plasma.level_number_density,
                plasma.time_explosion,
                plasma.stimulated_emission_factor,
            )

        opacity_state = OpacityState.from_legacy_plasma(plasma, tau_sobolev)

        return opacity_state

    def solve(self, plasma) -> OpacityState:
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
        atomic_data = plasma.atomic_data

        if self.disable_line_scattering:
            tau_sobolev = pd.DataFrame(
                np.zeros(
                    (
                        atomic_data.lines.shape[0],  # number of lines
                        plasma.number_density.shape[1],  # number of shells
                    ),
                    dtype=np.float64,
                ),
                index=atomic_data.lines.index,
            )
        else:
            tau_sobolev = calculate_sobolev_line_opacity(
                atomic_data.lines,
                plasma.level_number_density,
                plasma.time_explosion,
                plasma.stimulated_emission_factor,
            )
            beta_sobolev = calculate_beta_sobolev(
                tau_sobolev, self.previous_beta_sobolev
            )

        opacity_state = OpacityState.from_plasma(
            plasma, tau_sobolev, beta_sobolev
        )

        return opacity_state
