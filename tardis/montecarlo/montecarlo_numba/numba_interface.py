from enum import IntEnum

from numba import float64, int64, boolean
from numba.experimental import jitclass
import numpy as np

from astropy import units as u
from tardis import constants as const

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)


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
        r_inner: numpy.ndarray
        r_outer: numpy.ndarray
        time_explosion: float
        """
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.time_explosion = time_explosion


numba_plasma_spec = [
    ("electron_density", float64[:]),
    ("line_list_nu", float64[:]),
    ("tau_sobolev", float64[:, :]),
    ("transition_probabilities", float64[:, :]),
    ("line2macro_level_upper", int64[:]),
    ("macro_block_references", int64[:]),
    ("transition_type", int64[:]),
    ("destination_level_id", int64[:]),
    ("transition_line_id", int64[:]),
]


@jitclass(numba_plasma_spec)
class NumbaPlasma(object):
    def __init__(
        self,
        electron_density,
        line_list_nu,
        tau_sobolev,
        transition_probabilities,
        line2macro_level_upper,
        macro_block_references,
        transition_type,
        destination_level_id,
        transition_line_id,
    ):
        """
        Plasma for the Numba code
        Parameters
        ----------

        electron_density: numpy.array
        line_list_nu: numpy.array
        tau_sobolev: numpy.array
        transition_probabilities: numpy.array
        line2macro_level_upper: numpy.array
        macro_block_references: numpy.array
        transition_type: numpy.array
        destination_level_id: numpy.array
        transition_line_id: numpy.array
        """

        self.electron_density = electron_density
        self.line_list_nu = line_list_nu
        self.tau_sobolev = tau_sobolev

        #### Macro Atom transition probabilities
        self.transition_probabilities = transition_probabilities
        self.line2macro_level_upper = line2macro_level_upper

        self.macro_block_references = macro_block_references
        self.transition_type = transition_type

        # Destination level is not needed and/or generated for downbranch
        self.destination_level_id = destination_level_id
        self.transition_line_id = transition_line_id


def numba_plasma_initialize(plasma, line_interaction_type):
    """
    Initialize the NumbaPlasma object and copy over the data over from TARDIS Plasma

    Parameters
    -----------

    plasma: tardis.plasma.BasePlasma
    line_interaction_type: enum
    """
    electron_densities = plasma.electron_densities.values
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
        macro_block_references = plasma.atomic_data.macro_atom_references[
            "block_references"
        ].values
        transition_type = plasma.atomic_data.macro_atom_data[
            "transition_type"
        ].values

        # Destination level is not needed and/or generated for downbranch
        destination_level_id = plasma.atomic_data.macro_atom_data[
            "destination_level_idx"
        ].values
        transition_line_id = plasma.atomic_data.macro_atom_data[
            "lines_idx"
        ].values

    return NumbaPlasma(
        electron_densities,
        line_list_nu,
        tau_sobolev,
        transition_probabilities,
        line2macro_level_upper,
        macro_block_references,
        transition_type,
        destination_level_id,
        transition_line_id,
    )


packet_collection_spec = [
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
        packets_input_nu,
        packets_input_mu,
        packets_input_energy,
        packets_output_nu,
        packets_output_energy,
    ):
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
    ("idx", int64),
    ("number_of_vpackets", int64),
    ("length", int64),
    ("last_interaction_in_nu", float64),
    ("last_interaction_type", int64),
    ("last_interaction_in_id", int64),
    ("last_interaction_out_id", int64),
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
        self.number_of_vpackets = number_of_vpackets
        self.last_interaction_in_nu = 0.0
        self.last_interaction_type = -1
        self.last_interaction_in_id = -1
        self.last_interaction_out_id = -1
        self.idx = 0
        self.rpacket_index = rpacket_index
        self.length = temporary_v_packet_bins

    def set_properties(
        self,
        nu,
        energy,
        last_interaction_in_nu,
        last_interaction_type,
        last_interaction_in_id,
        last_interaction_out_id,
    ):
        if self.idx >= self.length:
            temp_length = self.length * 2 + self.number_of_vpackets
            temp_nus = np.empty(temp_length, dtype=np.float64)
            temp_energies = np.empty(temp_length, dtype=np.float64)
            temp_nus[: self.length] = self.nus
            temp_energies[: self.length] = self.energies

            self.nus = temp_nus
            self.energies = temp_energies
            self.length = temp_length

        self.nus[self.idx] = nu
        self.energies[self.idx] = energy
        self.last_interaction_type = last_interaction_type
        self.last_interaction_in_nu = last_interaction_in_nu
        self.last_interaction_in_id = last_interaction_in_id
        self.last_interaction_out_id = last_interaction_out_id
        self.idx += 1


estimators_spec = [
    ("j_estimator", float64[:]),
    ("nu_bar_estimator", float64[:]),
    ("j_blue_estimator", float64[:, :]),
    ("Edotlu_estimator", float64[:, :]),
]


@jitclass(estimators_spec)
class Estimators(object):
    def __init__(
        self, j_estimator, nu_bar_estimator, j_blue_estimator, Edotlu_estimator
    ):
        self.j_estimator = j_estimator
        self.nu_bar_estimator = nu_bar_estimator
        self.j_blue_estimator = j_blue_estimator
        self.Edotlu_estimator = Edotlu_estimator


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
