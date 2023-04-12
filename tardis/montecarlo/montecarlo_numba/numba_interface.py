from enum import IntEnum

from numba import float64, int64
from numba.experimental import jitclass
import numpy as np

from astropy import units as u
from tardis import constants as const

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)

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


numba_plasma_spec = [
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


@jitclass(numba_plasma_spec)
class NumbaPlasma(object):
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


def numba_plasma_initialize(plasma, line_interaction_type):
    """
    Initialize the NumbaPlasma object and copy over the data over from TARDIS Plasma

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
    if montecarlo_configuration.disable_line_scattering:
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

        if montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED:
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
    if montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED:
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
        ff_opacity_factor = plasma.ff_cooling_factor / np.sqrt(t_electrons)
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

    return NumbaPlasma(
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


packet_collection_spec = [
    ("packets_input_radius", float64[:]),
    ("packets_input_nu", float64[:]),
    ("packets_input_mu", float64[:]),
    ("packets_input_energy", float64[:]),
    ("packets_output_nu", float64[:]),
    ("packets_output_energy", float64[:]),
]


@jitclass(packet_collection_spec)
class PacketCollection(object):
    def __init__(
        self,
        packets_input_radius,
        packets_input_nu,
        packets_input_mu,
        packets_input_energy,
        packets_output_nu,
        packets_output_energy,
    ):
        self.packets_input_radius = packets_input_radius
        self.packets_input_nu = packets_input_nu
        self.packets_input_mu = packets_input_mu
        self.packets_input_energy = packets_input_energy
        self.packets_output_nu = packets_output_nu
        self.packets_output_energy = packets_output_energy


vpacket_collection_spec = [
    ("rpacket_index", int64),
    ("spectrum_frequency", float64[:]),
    ("v_packet_spawn_start_frequency", float64),
    ("v_packet_spawn_end_frequency", float64),
    ("nus", float64[:]),
    ("energies", float64[:]),
    ("initial_mus", float64[:]),
    ("initial_rs", float64[:]),
    ("idx", int64),
    ("number_of_vpackets", int64),
    ("length", int64),
    ("last_interaction_in_nu", float64[:]),
    ("last_interaction_type", int64[:]),
    ("last_interaction_in_id", int64[:]),
    ("last_interaction_out_id", int64[:]),
]


@jitclass(vpacket_collection_spec)
class VPacketCollection(object):
    def __init__(
        self,
        rpacket_index,
        spectrum_frequency,
        v_packet_spawn_start_frequency,
        v_packet_spawn_end_frequency,
        number_of_vpackets,
        temporary_v_packet_bins,
    ):
        self.spectrum_frequency = spectrum_frequency
        self.v_packet_spawn_start_frequency = v_packet_spawn_start_frequency
        self.v_packet_spawn_end_frequency = v_packet_spawn_end_frequency
        self.nus = np.empty(temporary_v_packet_bins, dtype=np.float64)
        self.energies = np.empty(temporary_v_packet_bins, dtype=np.float64)
        self.initial_mus = np.empty(temporary_v_packet_bins, dtype=np.float64)
        self.initial_rs = np.empty(temporary_v_packet_bins, dtype=np.float64)
        self.number_of_vpackets = number_of_vpackets
        self.last_interaction_in_nu = np.zeros(
            temporary_v_packet_bins, dtype=np.float64
        )
        self.last_interaction_type = -1 * np.ones(
            temporary_v_packet_bins, dtype=np.int64
        )
        self.last_interaction_in_id = -1 * np.ones(
            temporary_v_packet_bins, dtype=np.int64
        )
        self.last_interaction_out_id = -1 * np.ones(
            temporary_v_packet_bins, dtype=np.int64
        )
        self.idx = 0
        self.rpacket_index = rpacket_index
        self.length = temporary_v_packet_bins

    def set_properties(
        self,
        nu,
        energy,
        initial_mu,
        initial_r,
        last_interaction_in_nu,
        last_interaction_type,
        last_interaction_in_id,
        last_interaction_out_id,
    ):
        if self.idx >= self.length:
            temp_length = self.length * 2 + self.number_of_vpackets
            temp_nus = np.empty(temp_length, dtype=np.float64)
            temp_energies = np.empty(temp_length, dtype=np.float64)
            temp_initial_mus = np.empty(temp_length, dtype=np.float64)
            temp_initial_rs = np.empty(temp_length, dtype=np.float64)
            temp_last_interaction_in_nu = np.empty(
                temp_length, dtype=np.float64
            )
            temp_last_interaction_type = np.empty(temp_length, dtype=np.int64)
            temp_last_interaction_in_id = np.empty(temp_length, dtype=np.int64)
            temp_last_interaction_out_id = np.empty(temp_length, dtype=np.int64)

            temp_nus[: self.length] = self.nus
            temp_energies[: self.length] = self.energies
            temp_initial_mus[: self.length] = self.initial_mus
            temp_initial_rs[: self.length] = self.initial_rs
            temp_last_interaction_in_nu[
                : self.length
            ] = self.last_interaction_in_nu
            temp_last_interaction_type[
                : self.length
            ] = self.last_interaction_type
            temp_last_interaction_in_id[
                : self.length
            ] = self.last_interaction_in_id
            temp_last_interaction_out_id[
                : self.length
            ] = self.last_interaction_out_id

            self.nus = temp_nus
            self.energies = temp_energies
            self.initial_mus = temp_initial_mus
            self.initial_rs = temp_initial_rs
            self.last_interaction_in_nu = temp_last_interaction_in_nu
            self.last_interaction_type = temp_last_interaction_type
            self.last_interaction_in_id = temp_last_interaction_in_id
            self.last_interaction_out_id = temp_last_interaction_out_id
            self.length = temp_length

        self.nus[self.idx] = nu
        self.energies[self.idx] = energy
        self.initial_mus[self.idx] = initial_mu
        self.initial_rs[self.idx] = initial_r
        self.last_interaction_in_nu[self.idx] = last_interaction_in_nu
        self.last_interaction_type[self.idx] = last_interaction_type
        self.last_interaction_in_id[self.idx] = last_interaction_in_id
        self.last_interaction_out_id[self.idx] = last_interaction_out_id
        self.idx += 1


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
    ("interact_id", int64),
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
        interact_id : int
            Internal counter for the interactions that a particular RPacket undergoes
    """

    def __init__(self):
        self.length = montecarlo_configuration.INITIAL_TRACKING_ARRAY_LENGTH
        self.seed = np.int64(0)
        self.index = np.int64(0)
        self.status = np.empty(self.length, dtype=np.int64)
        self.r = np.empty(self.length, dtype=np.float64)
        self.nu = np.empty(self.length, dtype=np.float64)
        self.mu = np.empty(self.length, dtype=np.float64)
        self.energy = np.empty(self.length, dtype=np.float64)
        self.shell_id = np.empty(self.length, dtype=np.int64)
        self.interaction_type = np.empty(self.length, dtype=np.int64)
        self.interact_id = 0

    def track(self, r_packet):
        if self.interact_id >= self.length:
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
        self.status[self.interact_id] = r_packet.status
        self.r[self.interact_id] = r_packet.r
        self.nu[self.interact_id] = r_packet.nu
        self.mu[self.interact_id] = r_packet.mu
        self.energy[self.interact_id] = r_packet.energy
        self.shell_id[self.interact_id] = r_packet.current_shell_id
        self.interaction_type[self.interact_id] = r_packet.last_interaction_type
        self.interact_id += 1

    def finalize_array(self):
        self.status = self.status[: self.interact_id]
        self.r = self.r[: self.interact_id]
        self.nu = self.nu[: self.interact_id]
        self.mu = self.mu[: self.interact_id]
        self.energy = self.energy[: self.interact_id]
        self.shell_id = self.shell_id[: self.interact_id]
        self.interaction_type = self.interaction_type[: self.interact_id]


base_estimators_spec = [
    ("j_estimator", float64[:]),
    ("nu_bar_estimator", float64[:]),
    ("j_blue_estimator", float64[:, :]),
    ("Edotlu_estimator", float64[:, :]),
]

continuum_estimators_spec = [
    ("photo_ion_estimator", float64[:, :]),
    ("stim_recomb_estimator", float64[:, :]),
    ("bf_heating_estimator", float64[:, :]),
    ("stim_recomb_cooling_estimator", float64[:, :]),
    ("photo_ion_estimator_statistics", int64[:, :]),
]


@jitclass(base_estimators_spec + continuum_estimators_spec)
class Estimators(object):
    def __init__(
        self,
        j_estimator,
        nu_bar_estimator,
        j_blue_estimator,
        Edotlu_estimator,
        photo_ion_estimator,
        stim_recomb_estimator,
        bf_heating_estimator,
        stim_recomb_cooling_estimator,
        photo_ion_estimator_statistics,
    ):
        self.j_estimator = j_estimator
        self.nu_bar_estimator = nu_bar_estimator
        self.j_blue_estimator = j_blue_estimator
        self.Edotlu_estimator = Edotlu_estimator
        self.photo_ion_estimator = photo_ion_estimator
        self.stim_recomb_estimator = stim_recomb_estimator
        self.bf_heating_estimator = bf_heating_estimator
        self.stim_recomb_cooling_estimator = stim_recomb_cooling_estimator
        self.photo_ion_estimator_statistics = photo_ion_estimator_statistics

    def increment(self, other):

        self.j_estimator += other.j_estimator
        self.nu_bar_estimator += other.nu_bar_estimator
        self.j_blue_estimator += other.j_blue_estimator
        self.Edotlu_estimator += other.Edotlu_estimator
        self.photo_ion_estimator += other.photo_ion_estimator
        self.stim_recomb_estimator += other.stim_recomb_estimator
        self.bf_heating_estimator += other.bf_heating_estimator
        self.stim_recomb_cooling_estimator += (
            other.stim_recomb_cooling_estimator
        )
        self.photo_ion_estimator_statistics += (
            other.photo_ion_estimator_statistics
        )


def configuration_initialize(runner, number_of_vpackets):
    if runner.line_interaction_type == "macroatom":
        montecarlo_configuration.line_interaction_type = (
            LineInteractionType.MACROATOM
        )
    elif runner.line_interaction_type == "downbranch":
        montecarlo_configuration.line_interaction_type = (
            LineInteractionType.DOWNBRANCH
        )
    elif runner.line_interaction_type == "scatter":
        montecarlo_configuration.line_interaction_type = (
            LineInteractionType.SCATTER
        )
    else:
        raise ValueError(
            f'Line interaction type must be one of "macroatom",'
            f'"downbranch", or "scatter" but is '
            f"{runner.line_interaction_type}"
        )
    montecarlo_configuration.number_of_vpackets = number_of_vpackets
    montecarlo_configuration.temporary_v_packet_bins = number_of_vpackets
    montecarlo_configuration.full_relativity = runner.enable_full_relativity
    montecarlo_configuration.montecarlo_seed = runner.seed
    montecarlo_configuration.v_packet_spawn_start_frequency = (
        runner.virtual_spectrum_spawn_range.end.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
    montecarlo_configuration.v_packet_spawn_end_frequency = (
        runner.virtual_spectrum_spawn_range.start.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
    montecarlo_configuration.VPACKET_LOGGING = runner.virt_logging


# class TrackRPacket(object):
class LineInteractionType(IntEnum):
    SCATTER = 0
    DOWNBRANCH = 1
    MACROATOM = 2
