import numpy as np
import pandas as pd

from tardis.opacities.opacity_state import (
    OpacityState,
)
from tardis.transport.montecarlo.modes.nonhomologous.tau_sobolev import (
    calculate_beta_sobolev,
    calculate_beta_sobolev_directional,
    calculate_sobolev_line_opacity,
    calculate_sobolev_optical_depth_coefficient,
)


class OpacitySolver:
    velocity_gradient: np.ndarray
    line_interaction_type: str = "scatter"
    disable_line_scattering: bool = False

    def __init__(
        self,
        velocity_gradient,
        line_interaction_type="scatter",
        disable_line_scattering=False,
        velocity_over_radius=None,
    ):
        """Solver class for opacities

        Parameters
        ----------
        velocity_gradient : np.ndarray
        line_interaction_type: str
            "scatter", "downbranch", or "macroatom"
        disable_line_scattering: bool
        """
        self.velocity_gradient = velocity_gradient
        self.line_interaction_type = line_interaction_type
        self.disable_line_scattering = disable_line_scattering
        self.velocity_over_radius = velocity_over_radius

    def legacy_solve(self, plasma) -> OpacityState:
        """
        Solves the opacity state

        Parameters
        ----------
        plasma : tardis.plasma.BasePlasma
            legacy base plasma

        Returns
        -------
        OpacityState
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
        else:
            tau_sobolev = calculate_sobolev_line_opacity(
                plasma.atomic_data.lines,
                plasma.level_number_density,
                self.velocity_gradient,
                plasma.stimulated_emission_factor,
            )

        sobolev_optical_depth_coefficient = (
            calculate_sobolev_optical_depth_coefficient(
                plasma.atomic_data.lines,
                plasma.level_number_density,
                plasma.stimulated_emission_factor,
            )
        )
        if self.disable_line_scattering:
            sobolev_optical_depth_coefficient.iloc[:, :] = 0.0

        opacity_state = OpacityState.from_legacy_plasma(
            plasma,
            tau_sobolev,
            sobolev_optical_depth_coefficient=(
                sobolev_optical_depth_coefficient
            ),
        )
        if self.velocity_over_radius is not None:
            opacity_state.beta_sobolev = calculate_beta_sobolev_directional(
                sobolev_optical_depth_coefficient,
                self.velocity_gradient,
                self.velocity_over_radius,
            )

        return opacity_state

    def solve(self, plasma) -> OpacityState:
        """
        Solves the opacity state

        Parameters
        ----------
        plasma : tardis.plasma.BasePlasma
            legacy base plasma

        Returns
        -------
        OpacityState
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
                np.ones_like(tau_sobolev),
                index=tau_sobolev.index,
                columns=tau_sobolev.columns,
            )
        else:
            tau_sobolev = calculate_sobolev_line_opacity(
                plasma.atomic_data.lines,
                plasma.level_number_density,
                self.velocity_gradient,
                plasma.stimulated_emission_factor,
            )

        sobolev_optical_depth_coefficient = (
            calculate_sobolev_optical_depth_coefficient(
                plasma.atomic_data.lines,
                plasma.level_number_density,
                plasma.stimulated_emission_factor,
            )
        )
        if self.disable_line_scattering:
            sobolev_optical_depth_coefficient.iloc[:, :] = 0.0

        if self.velocity_over_radius is None:
            beta_sobolev = calculate_beta_sobolev(tau_sobolev)
        else:
            beta_sobolev = calculate_beta_sobolev_directional(
                sobolev_optical_depth_coefficient,
                self.velocity_gradient,
                self.velocity_over_radius,
            )

        opacity_state = OpacityState.from_plasma(
            plasma,
            tau_sobolev,
            beta_sobolev,
            sobolev_optical_depth_coefficient=(
                sobolev_optical_depth_coefficient
            ),
        )

        return opacity_state
