from tardis.opacities.tau_sobolev import calculate_sobolev_line_opacity
from tardis.opacities.opacity_state import (
    OpacityState,
)
from tardis.opacities.macro_atom.macroatom_solver import MacroAtomSolver
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
        if (
            self.line_interaction_type != "scatter" # TODO: Fix this
        ):  # Need a switch to use the continuum solver
            self.macro_atom_solver = MacroAtomSolver()
        else:
            self.macro_atom_solver = None

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
                        legacy_plasma.abundance.shape[1],  # number of shells
                    ),
                    dtype=np.float64,
                ),
                index=legacy_plasma.atomic_data.lines.index,
            )
        else:
            tau_sobolev = calculate_sobolev_line_opacity(
                legacy_plasma.atomic_data.lines,
                legacy_plasma.level_number_density,
                legacy_plasma.time_explosion,
                legacy_plasma.stimulated_emission_factor,
            )

        if self.line_interaction_type != "scatter":
            macroatom_state = self.macro_atom_solver.solve(
                legacy_plasma,
                atomic_data,
                tau_sobolev,
                legacy_plasma.stimulated_emission_factor,
            )
        else:
            macroatom_state = None

        opacity_state = OpacityState.from_legacy_plasma(
            legacy_plasma, tau_sobolev, macroatom_state
        )

        return opacity_state
