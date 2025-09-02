
import numba as nb
import numpy as np
from numba.experimental import jitclass

from tardis import constants as const
from tardis.transport.montecarlo.configuration import montecarlo_globals

C_SPEED_OF_LIGHT = const.c.to("cm/s").value


@jitclass
class OpacityStateNumba:
    electron_density: nb.float64[:]  # type: ignore[misc]
    t_electrons: nb.float64[:]  # type: ignore[misc]
    line_list_nu: nb.float64[:]  # type: ignore[misc]
    tau_sobolev: nb.float64[:, :]  # type: ignore[misc]
    transition_probabilities: nb.float64[:, :]  # type: ignore[misc]
    line2macro_level_upper: nb.int64[:]  # type: ignore[misc]
    macro_block_references: nb.int64[:]  # type: ignore[misc]
    transition_type: nb.int64[:]  # type: ignore[misc]
    destination_level_id: nb.int64[:]  # type: ignore[misc]
    transition_line_id: nb.int64[:]  # type: ignore[misc]
    bf_threshold_list_nu: nb.float64[:]  # type: ignore[misc]
    p_fb_deactivation: nb.float64[:, :]  # type: ignore[misc]
    photo_ion_nu_threshold_mins: nb.float64[:]  # type: ignore[misc]
    photo_ion_nu_threshold_maxs: nb.float64[:]  # type: ignore[misc]
    photo_ion_block_references: nb.int64[:]  # type: ignore[misc]
    chi_bf: nb.float64[:, :]  # type: ignore[misc]
    x_sect: nb.float64[:]  # type: ignore[misc]
    phot_nus: nb.float64[:]  # type: ignore[misc]
    ff_opacity_factor: nb.float64[:]  # type: ignore[misc]
    emissivities: nb.float64[:, :]  # type: ignore[misc]
    photo_ion_activation_idx: nb.int64[:]  # type: ignore[misc]
    k_packet_idx: nb.int64  # type: ignore[misc]

    def __init__(
        self,
        electron_density: np.ndarray,
        t_electrons: np.ndarray,
        line_list_nu: np.ndarray,
        tau_sobolev: np.ndarray,
        transition_probabilities: np.ndarray,
        line2macro_level_upper: np.ndarray,
        macro_block_references: np.ndarray,
        transition_type: np.ndarray,
        destination_level_id: np.ndarray,
        transition_line_id: np.ndarray,
        bf_threshold_list_nu: np.ndarray,
        p_fb_deactivation: np.ndarray,
        photo_ion_nu_threshold_mins: np.ndarray,
        photo_ion_nu_threshold_maxs: np.ndarray,
        photo_ion_block_references: np.ndarray,
        chi_bf: np.ndarray,
        x_sect: np.ndarray,
        phot_nus: np.ndarray,
        ff_opacity_factor: np.ndarray,
        emissivities: np.ndarray,
        photo_ion_activation_idx: np.ndarray,
        k_packet_idx: int,
    ) -> None:
        """
        Initialize Numba-compatible opacity state for Monte Carlo transport.

        Parameters
        ----------
        electron_density : numpy.ndarray
            Electron density in each shell [cm^-3].
        t_electrons : numpy.ndarray
            Electron temperature in each shell [K].
        line_list_nu : numpy.ndarray
            Frequencies of spectral lines [Hz].
        tau_sobolev : numpy.ndarray
            Sobolev optical depths for line transitions.
        transition_probabilities : numpy.ndarray
            Probabilities for macro atom transitions.
        line2macro_level_upper : numpy.ndarray
            Mapping from line indices to macro atom upper levels.
        macro_block_references : numpy.ndarray
            Block references for macro atom data.
        transition_type : numpy.ndarray
            Type identifiers for transitions.
        destination_level_id : numpy.ndarray
            Destination level indices for transitions.
        transition_line_id : numpy.ndarray
            Line indices for transitions.
        bf_threshold_list_nu : numpy.ndarray
            Bound-free threshold frequencies [Hz].
        p_fb_deactivation : numpy.ndarray
            Free-bound deactivation probabilities.
        photo_ion_nu_threshold_mins : numpy.ndarray
            Minimum photoionization threshold frequencies [Hz].
        photo_ion_nu_threshold_maxs : numpy.ndarray
            Maximum photoionization threshold frequencies [Hz].
        photo_ion_block_references : numpy.ndarray
            Block references for photoionization data.
        chi_bf : numpy.ndarray
            Bound-free absorption coefficients.
        x_sect : numpy.ndarray
            Photoionization cross sections [cm^2].
        phot_nus : numpy.ndarray
            Photoionization frequencies [Hz].
        ff_opacity_factor : numpy.ndarray
            Free-free opacity factors.
        emissivities : numpy.ndarray
            Emission coefficients for bound-free transitions.
        photo_ion_activation_idx : numpy.ndarray
            Indices for photoionization activation.
        k_packet_idx : int
            Index for k-packet handling.
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

    def __getitem__(self, i: slice) -> "OpacityStateNumba":
        """Get a shell or slice of shells of the attributes of the opacity state.

        Parameters
        ----------
        i : slice
            Shell slice. Will fail if slice is int since class only supports array types.

        Returns
        -------
        OpacityStateNumba
            A shallow copy of the current instance with sliced data.
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


def opacity_state_numba_initialize(
    plasma,
    line_interaction_type: str,
    disable_line_scattering: bool,
) -> OpacityStateNumba:
    """
    Initialize the OpacityStateNumba object and copy data from TARDIS Plasma.

    Parameters
    ----------
    plasma : tardis.plasma.BasePlasma
        The plasma object containing atomic and opacity data.
    line_interaction_type : str
        Type of line interaction ("scatter" or "macroatom").
    disable_line_scattering : bool
        Whether to disable line scattering by setting tau_sobolev to zero.

    Returns
    -------
    OpacityStateNumba
        Initialized opacity state for Monte Carlo transport.
    """
    electron_densities = plasma.electron_densities.values
    t_electrons = plasma.t_electrons
    line_list_nu = plasma.atomic_data.lines.nu.values
    tau_sobolev = np.ascontiguousarray(
        plasma.tau_sobolevs.values.copy(), dtype=np.float64
    )
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
        k_packet_idx = plasma.k_packet_idx
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
        k_packet_idx = -1

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


