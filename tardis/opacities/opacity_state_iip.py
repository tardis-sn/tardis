import numpy as np
from numba import float64, int64

from tardis.opacities.continuum.continuum_state import ContinuumState
from tardis.opacities.macro_atom.macroatom_state import MacroAtomState
from tardis.opacities.opacity_state import OpacityState
from tardis.opacities.opacity_state_numba_iip import OpacityStateNumbaIIP
from tardis.transport.montecarlo.configuration import montecarlo_globals


opacity_state_iip_spec = [
    ("electron_density", float64[:]),
    ("t_electrons", float64[:]),
    ("line_list_nu", float64[:]),
    ("tau_sobolev", float64[:, :]),
    ("transition_probabilities", float64[:, :]),
    ("line2macro_level_upper", int64[:]),
    ("macro_block_references", int64[:]),
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
    ("absorbing_markov_probabilities", float64[:, :]),
]


class OpacityStateIIP(OpacityState):
    """
    IIP-specific opacity state that extends OpacityState with support for
    absorbing Markov chain probabilities (matrix B) for direct-jump macroatom interactions.
    """

    def __init__(
        self,
        electron_density,
        t_electrons,
        line_list_nu,
        tau_sobolev,
        beta_sobolev,
        continuum_state,
        absorbing_markov_probabilities=None,
    ):
        """
        IIP Opacity State in Python

        Parameters
        ----------
        electron_density : pd.DataFrame
            Electron density per shell
        t_electrons : numpy.ndarray
            Electron temperature per shell
        line_list_nu : pd.DataFrame
            Line frequencies
        tau_sobolev : pd.DataFrame
            Expansion Optical Depths
        beta_sobolev : pd.DataFrame
            Modified expansion Optical Depths
        continuum_state : tardis.opacities.continuum.continuum_state.ContinuumState
            Continuum opacity data
        absorbing_markov_probabilities : pd.DataFrame, optional
            Matrix B: Absorbing probabilities of the Markov-chain macro atom.
            Indexed by source_level_idx, destination_level_idx per shell.
            Probability of being absorbed in destination_level_idx when
            starting from source_level_idx.
        """
        super().__init__(
            electron_density,
            t_electrons,
            line_list_nu,
            tau_sobolev,
            beta_sobolev,
            continuum_state,
        )
        self.absorbing_markov_probabilities = absorbing_markov_probabilities

    @classmethod
    def from_legacy_plasma(cls, plasma, tau_sobolev):
        """
        Generates an OpacityStateIIP object from a tardis BasePlasma

        Parameters
        ----------
        plasma : tardis.plasma.BasePlasma
            legacy base plasma
        tau_sobolev : pd.DataFrame
            Expansion Optical Depths

        Returns
        -------
        OpacityStateIIP
        """
        if hasattr(plasma, "photo_ion_cross_sections"):
            continuum_state = ContinuumState.from_legacy_plasma(plasma)
        else:
            continuum_state = None

        # Extract absorbing Markov probabilities if available
        absorbing_markov_probabilities = None
        if hasattr(plasma, "B"):
            absorbing_markov_probabilities = plasma.B

        return cls(
            plasma.electron_densities,
            plasma.t_electrons,
            plasma.atomic_data.lines.nu,
            tau_sobolev,
            plasma.beta_sobolev,
            continuum_state,
            absorbing_markov_probabilities,
        )

    @classmethod
    def from_plasma(cls, plasma, tau_sobolev, beta_sobolev):
        """
        Generates an OpacityStateIIP object from a tardis BasePlasma

        Parameters
        ----------
        plasma : tardis.plasma.BasePlasma
            legacy base plasma
        tau_sobolev : pd.DataFrame
            Expansion Optical Depths
        beta_sobolev : pd.DataFrame
            Modified expansion Optical Depths

        Returns
        -------
        OpacityStateIIP
        """
        if hasattr(plasma, "photo_ion_cross_sections"):
            continuum_state = ContinuumState.from_legacy_plasma(plasma)
        else:
            continuum_state = None

        # Extract absorbing Markov probabilities if available
        absorbing_markov_probabilities = None
        if hasattr(plasma, "B"):
            absorbing_markov_probabilities = plasma.B

        return cls(
            plasma.electron_densities,
            plasma.t_electrons,
            plasma.atomic_data.lines.nu,
            tau_sobolev,
            beta_sobolev,
            continuum_state,
            absorbing_markov_probabilities,
        )

    def to_numba(
        self,
        macro_atom_state: MacroAtomState,
        line_interaction_type,
    ) -> OpacityStateNumbaIIP:
        """
        Initialize the OpacityStateNumbaIIP object and copy over the data from OpacityStateIIP class

        Parameters
        ----------
        macro_atom_state : tardis.opacities.macro_atom.macroatom_state.MacroAtomState
            Macro atom state containing transition data
        line_interaction_type : enum
            Type of line interaction (scatter or macroatom)

        Returns
        -------
        OpacityStateNumbaIIP
            IIP-specific numba opacity state with absorbing Markov probabilities
        """
        electron_densities = self.electron_density.values
        t_electrons = self.t_electrons
        line_list_nu = self.line_list_nu.values

        # NOTE: Disabled line scattering is handled by the opacitystate solver
        tau_sobolev = np.ascontiguousarray(self.tau_sobolev, dtype=np.float64)

        if line_interaction_type == "scatter":
            # to adhere to data types, we must have an array of minimum size 1
            array_size = 1
            transition_probabilities = np.zeros(
                (array_size, array_size), dtype=np.float64
            )  # to adhere to data types
            line2macro_level_upper = np.zeros(array_size, dtype=np.int64)
            macro_block_references = np.zeros(array_size, dtype=np.int64)
            transition_type = np.zeros(array_size, dtype=np.int64)
            destination_level_id = np.zeros(array_size, dtype=np.int64)
            transition_line_id = np.zeros(array_size, dtype=np.int64)
        else:
            transition_probabilities = np.ascontiguousarray(
                macro_atom_state.transition_probabilities.values.copy(),
                dtype=np.float64,
            )
            line2macro_level_upper = (
                macro_atom_state.line2macro_level_upper.values
            )
            # TODO: Fix setting of block references for non-continuum mode

            macro_block_references = np.asarray(
                macro_atom_state.macro_block_references
            )

            transition_type = (
                macro_atom_state.transition_metadata.transition_type.values
            )

            # Destination level is not needed and/or generated for downbranch
            destination_level_id = macro_atom_state.transition_metadata.destination_level_idx.values
            transition_line_id = (
                macro_atom_state.transition_metadata.transition_line_idx.values
            )

        if montecarlo_globals.CONTINUUM_PROCESSES_ENABLED:
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
                self.continuum_state.photo_ion_activation_idx.values
            )
            k_packet_idx = np.int64(self.continuum_state.k_packet_idx)
        else:
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

        # Extract absorbing Markov probabilities if available
        if self.absorbing_markov_probabilities is not None:
            absorbing_markov_probabilities = np.ascontiguousarray(
                self.absorbing_markov_probabilities.values.copy(),
                dtype=np.float64,
            )
        else:
            # Initialize with empty array if not available
            # This allows the code to work even without the B matrix
            absorbing_markov_probabilities = np.zeros((1, 1), dtype=np.float64)

        return OpacityStateNumbaIIP(
            electron_densities,
            t_electrons,
            line_list_nu,
            tau_sobolev,
            transition_probabilities,
            line2macro_level_upper,
            macro_block_references,
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
