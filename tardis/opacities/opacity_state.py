from typing import Self

import numpy as np
import numpy.typing as npt
import pandas as pd
from numba import float64, int64

from tardis.opacities.continuum.continuum_state import ContinuumState
from tardis.opacities.macro_atom.macroatom_state import MacroAtomState
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.opacities.opacity_state_numba_iip import OpacityStateNumbaIIP
from tardis.plasma.base import BasePlasma
from tardis.transport.montecarlo.configuration import montecarlo_globals

opacity_state_spec = [
    ("electron_density", float64[:]),
    ("t_electrons", float64[:]),
    ("line_list_nu", float64[:]),
    ("tau_sobolev", float64[:, :]),
    ("transition_probabilities", float64[:, :]),
    ("line2macro_level_upper", int64[:]),
    ("macro_block_edge_index", int64[:]),
    ("transition_type", int64[:]),
    ("destination_level_id", int64[:]),
    ("transition_line_id", int64[:]),
    ("bf_threshold_list_nu", float64[:]),
    ("p_fb_deactivation", float64[:, :]),
    ("photo_ion_nu_threshold_mins", float64[:]),
    ("photo_ion_nu_threshold_maxs", float64[:]),
    ("photo_ion_block_references", int64[:]),
    ("chi_bf", float64[:, :]),
    ("x_sect", float64[:]),
    ("phot_nus", float64[:]),
    ("ff_opacity_factor", float64[:]),
    ("emissivities", float64[:, :]),
    ("photo_ion_activation_idx", int64[:]),
    ("k_packet_idx", int64),
]


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
        continuum_state: ContinuumState | None,
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
        continuum_state : tardis.opacities.continuum.continuum_state.ContinuumState or None
            Continuum quantities needed when continuum interactions are enabled.
        """
        self.electron_density = electron_density
        self.t_electrons = t_electrons
        self.line_list_nu = line_list_nu

        self.tau_sobolev = tau_sobolev

        self.beta_sobolev = beta_sobolev

        # Continuum Opacity Data
        self.continuum_state = continuum_state

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
        if hasattr(plasma, "photo_ion_cross_sections"):
            continuum_state = ContinuumState.from_legacy_plasma(plasma)
        else:
            continuum_state = None

        return cls(
            plasma.electron_densities,
            plasma.t_electrons,
            plasma.atomic_data.lines.nu,
            tau_sobolev,
            plasma.beta_sobolev,
            continuum_state,
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
        if continuum_state is None and hasattr(
            plasma, "photo_ion_cross_sections"
        ):
            continuum_state = ContinuumState.from_legacy_plasma(plasma)
        elif continuum_state is None:
            continuum_state = None

        return cls(
            plasma.electron_densities,
            plasma.t_electrons,
            plasma.atomic_data.lines.nu,
            tau_sobolev,
            beta_sobolev,
            continuum_state,
        )

    def to_numba(
        self,
        macro_atom_state: MacroAtomState | None,
        line_interaction_type: str,
    ) -> OpacityStateNumba | OpacityStateNumbaIIP:
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
        tardis.opacities.opacity_state_numba.OpacityStateNumba or tardis.opacities.opacity_state_numba_iip.OpacityStateNumbaIIP
            Array-backed opacity and interaction data used by transport.
        """
        electron_densities = self.electron_density.values
        t_electrons = self.t_electrons
        line_list_nu = self.line_list_nu.values

        # NOTE: Disabled line scattering is handled by the opacitystate solver
        tau_sobolev = np.ascontiguousarray(self.tau_sobolev, dtype=np.float64)

        # initialize the continuum attributes needed for the opacity state
        bf_threshold_list_nu = np.zeros(0, dtype=np.float64)
        p_fb_deactivation = np.zeros((0, 0), dtype=np.float64)
        photo_ion_nu_threshold_mins = np.zeros(0, dtype=np.float64)
        photo_ion_nu_threshold_maxs = np.zeros(0, dtype=np.float64)
        photo_ion_block_references = np.zeros(0, dtype=np.int64)
        chi_bf = np.zeros((0, 0), dtype=np.float64)
        x_sect = np.zeros(0, dtype=np.float64)
        phot_nus = np.zeros(0, dtype=np.float64)
        ff_opacity_factor = np.zeros(0, dtype=np.float64)
        emissivities = np.zeros((0, 0), dtype=np.float64)
        photo_ion_activation_idx = np.zeros(0, dtype=np.int64)
        k_packet_idx = np.int64(-1)

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

        elif (
            montecarlo_globals.CONTINUUM_PROCESSES_ENABLED
        ):  # continuum settings
            transition_probabilities = np.ascontiguousarray(
                macro_atom_state.normalized_deactivating_probs.values.copy(),
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

            bf_threshold_list_nu = (
                self.continuum_state.bf_threshold_list_nu.values
            )
            p_fb_deactivation = np.ascontiguousarray(
                self.continuum_state.p_fb_deactivation.values.copy(),
                dtype=np.float64,
            )

            phot_nus = self.continuum_state.phot_nus
            photo_ion_block_references = (
                self.continuum_state.photo_ion_block_references
            )
            photo_ion_nu_threshold_mins = (
                self.continuum_state.photo_ion_nu_threshold_mins.values
            )
            photo_ion_nu_threshold_maxs = (
                self.continuum_state.photo_ion_nu_threshold_maxs.values
            )

            chi_bf = self.continuum_state.chi_bf.values
            x_sect = self.continuum_state.x_sect.values

            phot_nus = phot_nus.values
            ff_opacity_factor = (
                self.continuum_state.ff_cooling_factor / np.sqrt(t_electrons)
            ).astype(np.float64)
            emissivities = self.continuum_state.emissivities.values
            photo_ion_activation_idx = (
                self.continuum_state.photo_ion_activation_idx.to_numpy(
                    dtype=np.int64
                )
            )
            k_packet_idx = np.int64(macro_atom_state.k_packet_idx)
            absorbing_markov_probabilities = (
                macro_atom_state.absorbing_probability_matrix
            )
            return OpacityStateNumbaIIP(
                electron_densities,
                t_electrons,
                line_list_nu,
                tau_sobolev,
                transition_probabilities,
                line2macro_level_upper,
                macro_block_edge_index,
                transition_type,
                destination_level_id,
                transition_line_id,
                bf_threshold_list_nu,
                p_fb_deactivation,
                photo_ion_nu_threshold_mins,
                photo_ion_nu_threshold_maxs,
                photo_ion_block_references,
                chi_bf,
                x_sect,
                phot_nus,
                ff_opacity_factor,
                emissivities,
                photo_ion_activation_idx,
                k_packet_idx,
                absorbing_markov_probabilities,
            )
        else:
            # Not continuum
            transition_probabilities = np.ascontiguousarray(
                macro_atom_state.transition_probabilities.values.copy(),
                dtype=np.float64,
            )
            line2macro_level_upper = (
                macro_atom_state.line2macro_level_upper.values
            )
            # TODO: Fix setting of block references for non-continuum mode

            macro_block_edge_index = np.asarray(
                macro_atom_state.macro_block_edge_index
            )

            transition_type = (
                macro_atom_state.transition_metadata.transition_type.values
            )

            # Destination level is not needed and/or generated for downbranch
            destination_level_id = macro_atom_state.transition_metadata.destination_level_idx.values
            transition_line_id = (
                macro_atom_state.transition_metadata.transition_line_idx.values
            )

        return OpacityStateNumba(
            electron_densities,
            t_electrons,
            line_list_nu,
            tau_sobolev,
            transition_probabilities,
            line2macro_level_upper,
            macro_block_edge_index,
            transition_type,
            destination_level_id,
            transition_line_id,
            bf_threshold_list_nu,
            p_fb_deactivation,
            photo_ion_nu_threshold_mins,
            photo_ion_nu_threshold_maxs,
            photo_ion_block_references,
            chi_bf,
            x_sect,
            phot_nus,
            ff_opacity_factor,
            emissivities,
            photo_ion_activation_idx,
            k_packet_idx,
        )
