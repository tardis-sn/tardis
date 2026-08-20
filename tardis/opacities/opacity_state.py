from typing import Self

import numpy as np
import numpy.typing as npt
import pandas as pd

from tardis.opacities.macro_atom.macroatom_state import MacroAtomState
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.plasma.base import BasePlasma


class OpacityState:
    """Store Python-native line and continuum opacity data for one iteration.

    The state preserves labelled plasma data for the MC solver and formal-integral.
    Use :meth:`to_numba` to produce the Numba-compatible transport
    representation.
    """

    def __init__(
        self,
        electron_density: pd.DataFrame,
        t_electrons: npt.NDArray[np.float64],
        line_list_nu: pd.Series,
        tau_sobolev: pd.DataFrame,
        beta_sobolev: pd.DataFrame | None,
    ) -> None:
        """
        Initialize the Python-native opacity state.

        Parameters
        ----------
        electron_density : pd.DataFrame
        t_electrons : numpy.ndarray
            Electron temperatures in each shell [K].
        line_list_nu : pd.DataFrame
        tau_sobolev : pd.DataFrame
        beta_sobolev : pd.DataFrame or None
            Sobolev escape probabilities for each line and shell.
        """
        self.electron_density = electron_density
        self.t_electrons = t_electrons
        self.line_list_nu = line_list_nu

        self.tau_sobolev = tau_sobolev

        self.beta_sobolev = beta_sobolev

    @classmethod
    def from_legacy_plasma(
        cls,
        plasma: BasePlasma,
        tau_sobolev: pd.DataFrame,
    ) -> Self:
        """
        Construct an opacity state from a legacy plasma object.

        Parameters
        ----------
        plasma : tardis.plasma.BasePlasma
            Plasma object containing the line and continuum quantities.
        tau_sobolev : pd.DataFrame
            Sobolev optical depths for each line and shell.

        Returns
        -------
        OpacityState
            Python-native opacity state.
        """
        return cls(
            plasma.electron_densities,
            plasma.t_electrons,
            plasma.atomic_data.lines.nu,
            tau_sobolev,
            plasma.beta_sobolev,
        )

    @classmethod
    def from_plasma(
        cls,
        plasma: BasePlasma,
        tau_sobolev: pd.DataFrame,
        beta_sobolev: pd.DataFrame | None,
    ) -> Self:
        """
        Construct an opacity state from a plasma object.

        Parameters
        ----------
        plasma : tardis.plasma.base.BasePlasma
            Plasma object containing the line and continuum quantities.
        tau_sobolev : pd.DataFrame
            Sobolev optical depths for each line and shell.
        beta_sobolev : pd.DataFrame or None
            Sobolev escape probabilities for each line and shell.

        Returns
        -------
        OpacityState
            Python-native opacity state.
        """
        return cls(
            plasma.electron_densities,
            plasma.t_electrons,
            plasma.atomic_data.lines.nu,
            tau_sobolev,
            beta_sobolev,
        )

    def to_numba(
        self,
        macro_atom_state: MacroAtomState | None,
        line_interaction_type: str,
        continuum_processes_enabled: bool = False,
    ) -> OpacityStateNumba:
        """
        Convert this state to the Numba-compatible transport representation.

        Parameters
        ----------
        macro_atom_state : tardis.opacities.macro_atom.macroatom_state.MacroAtomState or None
            Macro-atom transition data. It is required unless
            ``line_interaction_type`` is ``"scatter"``.
        line_interaction_type : str
            Configured line-interaction mode.

        Returns
        -------
        tardis.opacities.opacity_state_numba.OpacityStateNumba
            Array-backed opacity and interaction data used by transport.
        """
        if line_interaction_type == "scatter":
            # to adhere to data types, we must have an array of minimum size 1
            array_size = 1
            transition_probabilities = np.zeros(
                (array_size, array_size), dtype=np.float64
            )  # to adhere to data types
            line2macro_level_upper = np.zeros(array_size, dtype=np.int64)
            macro_block_edge_index = np.zeros(array_size, dtype=np.int64)
            transition_type = np.zeros(array_size, dtype=np.int64)
            destination_level_id = np.zeros(array_size, dtype=np.int64)
            transition_line_id = np.zeros(array_size, dtype=np.int64)

        else:
            transition_probabilities = np.ascontiguousarray(
                (
                    macro_atom_state.normalized_deactivating_probs
                    if continuum_processes_enabled
                    else macro_atom_state.transition_probabilities
                ).values.copy(),
                dtype=np.float64,
            )
            line2macro_level_upper = (
                macro_atom_state.line2macro_level_upper.values
            )
            macro_block_edge_index = np.asarray(
                macro_atom_state.macro_block_edge_index
            )
            transition_type = (
                macro_atom_state.transition_metadata.transition_type.values
            )
            destination_level_id = macro_atom_state.transition_metadata.destination_level_idx.values
            transition_line_id = (
                macro_atom_state.transition_metadata.transition_line_idx.values
            )
        return OpacityStateNumba(
            self.electron_density.values,
            self.t_electrons,
            self.line_list_nu.values,
            np.ascontiguousarray(self.tau_sobolev, dtype=np.float64),
            transition_probabilities,
            line2macro_level_upper,
            macro_block_edge_index,
            transition_type,
            destination_level_id,
            transition_line_id,
        )
