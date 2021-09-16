from enum import IntEnum

from numba import float64, int64, boolean, njit
from numba.experimental import jitclass
import numpy as np

from astropy import units as u
from tardis import constants as const

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
from tardis.montecarlo.montecarlo_numba import njit_dict_no_parallel

C_SPEED_OF_LIGHT = const.c.to("cm/s").value

numba_model_spec = [
    ("r_inner", float64[:]),
    ("r_outer", float64[:]),
    ("time_explosion", float64),
]


@jitclass(numba_model_spec)
class NumbaModel(object):
    def __init__(self, r_inner, r_outer, time_explosion):
        """
        Model for the Numba mode

        Parameters
        ----------
        r_inner : numpy.ndarray
        r_outer : numpy.ndarray
        time_explosion : float
        """
        self.r_inner = r_inner
        self.r_outer = r_outer
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
        p_fb_deactivation
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
    p_fb_deactivation = np.ascontiguousarray(
            plasma.p_fb_deactivation.values.copy(), dtype=np.float64)
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
        # macro_block_references = plasma.atomic_data.macro_atom_references[
        #     "block_references"
        # ].values
        # TODO: Fix setting of block references for non-continuum mode
        macro_block_references = plasma.macro_block_references
        transition_type = plasma.macro_atom_data["transition_type"].values

        # Destination level is not needed and/or generated for downbranch
        destination_level_id = plasma.macro_atom_data[
            "destination_level_idx"
        ].values
        transition_line_id = plasma.macro_atom_data["lines_idx"].values
    if not plasma.continuum_interaction_species.empty:
        bf_threshold_list_nu = plasma.nu_i.loc[
            plasma.level2continuum_idx.index
        ].values
    else:
        bf_threshold_list_nu = np.zeros(0, dtype=np.int64)

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
        p_fb_deactivation
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



def create_continuum_class(plasma):
    """called before mainloop"""

    chi_continuum_calculator = plasma.chi_continuum_calculator
    nu_fb_sampler = plasma.nu_fb_sampler
    nu_ff_sampler = plasma.nu_ff_sampler
    get_macro_activation_idx = plasma.determine_continuum_macro_activation_idx
    continuum_spec = [
            ("chi_bf_tot", float64),
            ("chi_bf_contributions", float64[:]),
            ("current_continua", int64[:]),
            ("x_sect_bfs", float64[:]),
            ("chi_ff", float64),
    ]
    @jitclass(continuum_spec)
    class Continuum(object):

        def __init__(self):

            self.chi_bf_tot = 0.0
            self.chi_bf_contributions = np.empty(0, dtype=float64)
            self.current_continua = np.empty(0, dtype=int64)
            self.x_sect_bfs = np.empty(0, dtype=float64)
            self.chi_ff = 0.0

        def calculate(self, nu, shell):

            (
            self.chi_bf_tot,
            self.chi_bf_contributions,
            self.current_continua,
            self.x_sect_bfs,
            self.chi_ff,
            ) = chi_continuum_calculator(nu, shell)

        def sample_nu_free_bound(self, shell, continuum_id):

            return nu_fb_sampler(shell, continuum_id)

        def sample_nu_free_free(self, shell):

            return nu_ff_sampler(shell)

        def determine_macro_activation_idx(self, nu, shell):

            idx = get_macro_activation_idx(
                    nu, self.chi_bf_tot, self.chi_ff, 
                    self.chi_bf_contributions, self.current_continua
                    )
            return idx


    @njit(**njit_dict_no_parallel)
    def continuum_constructor():
        return Continuum()

    return continuum_constructor




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
    montecarlo_configuration.single_packet_seed = runner.single_packet_seed
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
