"""
IIP-specific opacity solver that creates OpacityStateIIP with absorbing Markov probabilities.
"""

import numpy as np
import pandas as pd

from tardis.opacities.opacity_solver import OpacitySolver
from tardis.opacities.opacity_state_iip import OpacityStateIIP
from tardis.opacities.tau_sobolev import (
    calculate_beta_sobolev,
    calculate_sobolev_line_opacity,
)


class OpacitySolverIIP(OpacitySolver):
    """
    IIP-specific opacity solver that extends OpacitySolver to create
    OpacityStateIIP instances with absorbing Markov chain probabilities.
    """

    def solve(self, plasma) -> OpacityStateIIP:
        """
        Solves the opacity state for IIP mode, including absorbing Markov probabilities.

        Parameters
        ----------
        plasma : tardis.plasma.BasePlasma or tardis.iip_plasma.base.LegacyPlasmaArray
            Plasma object from which to extract opacity data

        Returns
        -------
        OpacityStateIIP
            IIP opacity state with absorbing Markov chain probabilities if available
        """
        if self.disable_line_scattering:
            tau_sobolev = pd.DataFrame(
                np.zeros(
                    (
                        plasma.atomic_data.lines.shape[0],  # number of lines
                        plasma.number_density.shape[1],  # number of shells
                    ),
                    dtype=np.float64,
                ),
                index=plasma.atomic_data.lines.index,
            )
            beta_sobolev = pd.DataFrame(
                np.zeros_like(tau_sobolev.values),
                index=tau_sobolev.index,
                columns=tau_sobolev.columns,
            )
        else:
            tau_sobolev = calculate_sobolev_line_opacity(
                plasma.atomic_data.lines,
                plasma.level_number_density,
                plasma.time_explosion,
                plasma.stimulated_emission_factor,
            )
            beta_sobolev = calculate_beta_sobolev(tau_sobolev)

        opacity_state = OpacityStateIIP.from_plasma(
            plasma, tau_sobolev, beta_sobolev
        )

        return opacity_state
