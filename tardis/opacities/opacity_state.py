import numpy as np
from numba import float64, int64
from numba.experimental import jitclass

from tardis.opacities.tau_sobolev import calculate_sobolev_line_opacity
from tardis.transport.montecarlo.configuration import montecarlo_globals
from tardis.opacities.continuum.continuum_state import ContinuumState
from tardis.opacities.macro_atom.macroatom_state import MacroAtomState


class OpacityState:
    def __init__(
        self,
        electron_density,
        t_electrons,
        line_list_nu,
        tau_sobolev,
        continuum_state,
    ):
        """
        Opacity State in Python

        Parameters
        ----------
        electron_density : pd.DataFrame
        t_electrons : numpy.ndarray
        line_list_nu : pd.DataFrame
        tau_sobolev : pd.DataFrame
        continuum_state: tardis.opacities.continuum.continuum_state.ContinuumState
        """
        self.electron_density = electron_density
        self.t_electrons = t_electrons
        self.line_list_nu = line_list_nu

        self.tau_sobolev = tau_sobolev

        # Continuum Opacity Data
        self.continuum_state = continuum_state

    @classmethod
    def from_legacy_plasma(cls, plasma, tau_sobolev):
        """
        Generates an OpacityStatePython object from a tardis BasePlasma

        Parameters
        ----------
        plasma : tarids.plasma.BasePlasma
            legacy base plasma
        tau_sobolev : pd.DataFrame
            Expansion Optical Depths

        Returns
        -------
        OpacityStatePython
        """

        if hasattr(plasma, "photo_ion_cross_sections"):
            continuum_state = ContinuumState.from_legacy_plasma(plasma)
        else:
            continuum_state = None

        atomic_data = plasma.atomic_data

        return cls(
            plasma.electron_densities,
            plasma.t_electrons,
            atomic_data.lines.nu,
            tau_sobolev,
            continuum_state,
        )


opacity_state_spec = [
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
]


@jitclass(opacity_state_spec)
class OpacityStateNumba:
    def __init__(
        self,
        electron_density,
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
    ):
        """
        Plasma for the Numba code

        Parameters
        ----------
        electron_density : numpy.ndarray
        t_electrons : numpy.ndarray
        line_list_nu : numpy.ndarray
        tau_sobolev : numpy.ndarray
        transition_probabilities : numpy.ndarray
        line2macro_level_upper : numpy.ndarray
        macro_block_references : numpy.ndarray
        transition_type : numpy.ndarray
        destination_level_id : numpy.ndarray
        transition_line_id : numpy.ndarray
        bf_threshold_list_nu : numpy.ndarray
        """
        self.electron_density = electron_density
        self.t_electrons = t_electrons
        self.line_list_nu = line_list_nu
        self.tau_sobolev = tau_sobolev
        self.bf_threshold_list_nu = bf_threshold_list_nu

        #### Macro Atom transition probabilities
        self.transition_probabilities = transition_probabilities
        self.line2macro_level_upper = line2macro_level_upper

        self.macro_block_references = macro_block_references
        self.transition_type = transition_type

        # Destination level is not needed and/or generated for downbranch
        self.destination_level_id = destination_level_id
        self.transition_line_id = transition_line_id
        self.p_fb_deactivation = p_fb_deactivation

        # Continuum Opacity Data
        self.photo_ion_nu_threshold_mins = photo_ion_nu_threshold_mins
        self.photo_ion_nu_threshold_maxs = photo_ion_nu_threshold_maxs

        self.photo_ion_block_references = photo_ion_block_references
        self.chi_bf = chi_bf
        self.x_sect = x_sect
        self.phot_nus = phot_nus
        self.ff_opacity_factor = ff_opacity_factor
        self.emissivities = emissivities
        self.photo_ion_activation_idx = photo_ion_activation_idx
        self.k_packet_idx = k_packet_idx

    def __getitem__(self, i: slice):
        """Get a shell or slice of shells of the attributes of the opacity state

        Args:
            i (slice): shell slice.  Will fail if slice is int since class only supports array types

        Returns:
            OpacityState : a shallow copy of the current instance
        """
        # NOTE: This currently will not work with continuum processes since it does not slice those arrays
        return OpacityStateNumba(
            self.electron_density[i],
            self.t_electrons[i],
            self.line_list_nu,
            self.tau_sobolev[:, i],
            self.transition_probabilities[:, i],
            self.line2macro_level_upper,
            self.macro_block_references,
            self.transition_type,
            self.destination_level_id,
            self.transition_line_id,
            self.bf_threshold_list_nu,
            self.p_fb_deactivation,
            self.photo_ion_nu_threshold_mins,
            self.photo_ion_nu_threshold_maxs,
            self.photo_ion_block_references,
            self.chi_bf,
            self.x_sect,
            self.phot_nus,
            self.ff_opacity_factor,
            self.emissivities,
            self.photo_ion_activation_idx,
            self.k_packet_idx,
        )


def opacity_state_to_numba(
    opacity_state: OpacityState,
    macro_atom_state: MacroAtomState,
    line_interaction_type,
) -> OpacityStateNumba:
    """
    Initialize the OpacityStateNumba object and copy over the data over from OpacityState class

    Parameters
    ----------
    opacity_state : tardis.opacities.opacity_state.OpacityState
    line_interaction_type : enum
    """

    electron_densities = opacity_state.electron_density.values
    t_electrons = opacity_state.t_electrons
    line_list_nu = opacity_state.line_list_nu.values

    # NOTE: Disabled line scattering is handled by the opacitystate solver
    tau_sobolev = np.ascontiguousarray(
        opacity_state.tau_sobolev, dtype=np.float64
    )

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
        line2macro_level_upper = macro_atom_state.line2macro_level_upper
        # TODO: Fix setting of block references for non-continuum mode

        macro_block_references = np.asarray(
            macro_atom_state.macro_block_references
        )

        transition_type = macro_atom_state.transition_type.values

        # Destination level is not needed and/or generated for downbranch
        destination_level_id = macro_atom_state.destination_level_id.values
        transition_line_id = macro_atom_state.transition_line_id.values

    if montecarlo_globals.CONTINUUM_PROCESSES_ENABLED:
        bf_threshold_list_nu = (
            opacity_state.continuum_state.bf_threshold_list_nu.values
        )
        p_fb_deactivation = np.ascontiguousarray(
            opacity_state.continuum_state.p_fb_deactivation.values.copy(),
            dtype=np.float64,
        )

        phot_nus = opacity_state.continuum_state.phot_nus
        photo_ion_block_references = (
            opacity_state.continuum_state.photo_ion_block_references
        )
        photo_ion_nu_threshold_mins = (
            opacity_state.continuum_state.photo_ion_nu_threshold_mins.values
        )
        photo_ion_nu_threshold_maxs = (
            opacity_state.continuum_state.photo_ion_nu_threshold_maxs.values
        )

        chi_bf = opacity_state.continuum_state.chi_bf.values
        x_sect = opacity_state.continuum_state.x_sect.values

        phot_nus = phot_nus.values
        ff_opacity_factor = (
            opacity_state.continuum_state.ff_cooling_factor
            / np.sqrt(t_electrons)
        ).astype(np.float64)
        emissivities = opacity_state.continuum_state.emissivities.values
        photo_ion_activation_idx = (
            opacity_state.continuum_state.photo_ion_activation_idx.values
        )
        k_packet_idx = np.int64(opacity_state.continuum_state.k_packet_idx)
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

    return OpacityStateNumba(
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
    )


def opacity_state_initialize(
    plasma,
    line_interaction_type,
    disable_line_scattering,
):
    """
    Initialize the OpacityState object and copy over the data over from TARDIS Plasma

    Parameters
    ----------
    plasma : tardis.plasma.BasePlasma
    line_interaction_type : enum
    """
    electron_densities = plasma.electron_densities.values
    t_electrons = plasma.t_electrons
    line_list_nu = plasma.atomic_data.lines.nu.values

    tau_sobolev_df = calculate_sobolev_line_opacity(
        plasma.atomic_data.lines,
        plasma.level_number_density,
        plasma.time_explosion,
        plasma.stimulated_emission_factor,
    )

    tau_sobolev = np.ascontiguousarray(tau_sobolev_df, dtype=np.float64)

    if disable_line_scattering:
        tau_sobolev *= 0

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
            plasma.transition_probabilities.values.copy(), dtype=np.float64
        )
        line2macro_level_upper = (
            plasma.atomic_data.lines_upper2macro_reference_idx
        )
        # TODO: Fix setting of block references for non-continuum mode

        if montecarlo_globals.CONTINUUM_PROCESSES_ENABLED:
            macro_block_references = plasma.macro_block_references
        else:
            macro_block_references = plasma.atomic_data.macro_atom_references[
                "block_references"
            ].values
        transition_type = plasma.macro_atom_data["transition_type"].values

        # Destination level is not needed and/or generated for downbranch
        destination_level_id = plasma.macro_atom_data[
            "destination_level_idx"
        ].values
        transition_line_id = plasma.macro_atom_data["lines_idx"].values
    if montecarlo_globals.CONTINUUM_PROCESSES_ENABLED:
        bf_threshold_list_nu = plasma.nu_i.loc[
            plasma.level2continuum_idx.index
        ].values
        p_fb_deactivation = np.ascontiguousarray(
            plasma.p_fb_deactivation.values.copy(), dtype=np.float64
        )

        phot_nus = plasma.photo_ion_cross_sections.nu.loc[
            plasma.level2continuum_idx.index
        ]
        photo_ion_block_references = np.pad(
            phot_nus.groupby(level=[0, 1, 2], sort=False)
            .count()
            .values.cumsum(),
            [1, 0],
        )
        photo_ion_nu_threshold_mins = (
            phot_nus.groupby(level=[0, 1, 2], sort=False).first().values
        )
        photo_ion_nu_threshold_maxs = (
            phot_nus.groupby(level=[0, 1, 2], sort=False).last().values
        )

        chi_bf = plasma.chi_bf.loc[plasma.level2continuum_idx.index].values
        x_sect = plasma.photo_ion_cross_sections.x_sect.loc[
            plasma.level2continuum_idx.index
        ].values

        phot_nus = phot_nus.values
        ff_opacity_factor = (
            plasma.ff_cooling_factor / np.sqrt(t_electrons)
        ).astype(np.float64)
        emissivities = plasma.fb_emission_cdf.loc[
            plasma.level2continuum_idx.index
        ].values
        photo_ion_activation_idx = plasma.photo_ion_idx.loc[
            plasma.level2continuum_idx.index, "destination_level_idx"
        ].values
        k_packet_idx = np.int64(plasma.k_packet_idx)
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

    return OpacityStateNumba(
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
    )
