from enum import IntEnum

from numba import float64, int64
from numba.experimental import jitclass
import numpy as np

from tardis import constants as const


C_SPEED_OF_LIGHT = const.c.to("cm/s").value


numba_model_spec = [
    ("time_explosion", float64),
]


@jitclass(numba_model_spec)
class NumbaModel(object):
    def __init__(self, time_explosion):
        """
        Model for the Numba mode

        Parameters
        ----------
        time_explosion : float
        """
        self.time_explosion = time_explosion


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
class OpacityState(object):
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


def opacity_state_initialize(
    plasma,
    line_interaction_type,
    disable_line_scattering,
    continuum_processes_enabled,
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

        if continuum_processes_enabled:
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
    if continuum_processes_enabled:
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

    return OpacityState(
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


rpacket_tracker_spec = [
    ("length", int64),
    ("seed", int64),
    ("index", int64),
    ("status", int64[:]),
    ("r", float64[:]),
    ("nu", float64[:]),
    ("mu", float64[:]),
    ("energy", float64[:]),
    ("shell_id", int64[:]),
    ("interaction_type", int64[:]),
    ("num_interactions", int64),
]


@jitclass(rpacket_tracker_spec)
class RPacketTracker(object):
    """
    Numba JITCLASS for storing the information for each interaction a RPacket instance undergoes.

    Parameters
    ----------
        length : int
            Length of the initial array that is instantiated
        seed : int
            Seed for each RPacket
        index : int
            Index position of each RPacket
        status : int
            Current status of the RPacket as per interactions
        r : float
            Radius of the shell where the RPacket is present
        nu : float
            Luminosity of the RPacket
        mu : float
            Cosine of the angle made by the direction of movement of the RPacket from its original direction
        energy : float
            Energy possessed by the RPacket at a particular shell
        shell_id : int
            Current Shell No in which the RPacket is present
        interaction_type: int
            Type of interaction the rpacket undergoes
        num_interactions : int
            Internal counter for the interactions that a particular RPacket undergoes
    """

    def __init__(self, length):
        self.length = length
        self.seed = np.int64(0)
        self.index = np.int64(0)
        self.status = np.empty(self.length, dtype=np.int64)
        self.r = np.empty(self.length, dtype=np.float64)
        self.nu = np.empty(self.length, dtype=np.float64)
        self.mu = np.empty(self.length, dtype=np.float64)
        self.energy = np.empty(self.length, dtype=np.float64)
        self.shell_id = np.empty(self.length, dtype=np.int64)
        self.interaction_type = np.empty(self.length, dtype=np.int64)
        self.num_interactions = 0

    def track(self, r_packet):
        if self.num_interactions >= self.length:
            temp_length = self.length * 2
            temp_status = np.empty(temp_length, dtype=np.int64)
            temp_r = np.empty(temp_length, dtype=np.float64)
            temp_nu = np.empty(temp_length, dtype=np.float64)
            temp_mu = np.empty(temp_length, dtype=np.float64)
            temp_energy = np.empty(temp_length, dtype=np.float64)
            temp_shell_id = np.empty(temp_length, dtype=np.int64)
            temp_interaction_type = np.empty(temp_length, dtype=np.int64)

            temp_status[: self.length] = self.status
            temp_r[: self.length] = self.r
            temp_nu[: self.length] = self.nu
            temp_mu[: self.length] = self.mu
            temp_energy[: self.length] = self.energy
            temp_shell_id[: self.length] = self.shell_id
            temp_interaction_type[: self.length] = self.interaction_type

            self.status = temp_status
            self.r = temp_r
            self.nu = temp_nu
            self.mu = temp_mu
            self.energy = temp_energy
            self.shell_id = temp_shell_id
            self.interaction_type = temp_interaction_type
            self.length = temp_length

        self.index = r_packet.index
        self.seed = r_packet.seed
        self.status[self.num_interactions] = r_packet.status
        self.r[self.num_interactions] = r_packet.r
        self.nu[self.num_interactions] = r_packet.nu
        self.mu[self.num_interactions] = r_packet.mu
        self.energy[self.num_interactions] = r_packet.energy
        self.shell_id[self.num_interactions] = r_packet.current_shell_id
        self.interaction_type[
            self.num_interactions
        ] = r_packet.last_interaction_type
        self.num_interactions += 1

    def finalize_array(self):
        self.status = self.status[: self.num_interactions]
        self.r = self.r[: self.num_interactions]
        self.nu = self.nu[: self.num_interactions]
        self.mu = self.mu[: self.num_interactions]
        self.energy = self.energy[: self.num_interactions]
        self.shell_id = self.shell_id[: self.num_interactions]
        self.interaction_type = self.interaction_type[: self.num_interactions]


class LineInteractionType(IntEnum):
    SCATTER = 0
    DOWNBRANCH = 1
    MACROATOM = 2
