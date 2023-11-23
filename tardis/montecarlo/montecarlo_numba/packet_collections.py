import numpy as np
from numba.experimental import jitclass
from numba import float64, int64

packet_collection_spec = [
    ("initial_radii", float64[:]),
    ("initial_nus", float64[:]),
    ("initial_mus", float64[:]),
    ("initial_energies", float64[:]),
    ("time_of_simulation", float64),
    ("radiation_field_luminosity", float64),  #
    ("output_nus", float64[:]),
    ("output_energies", float64[:]),
]


@jitclass(packet_collection_spec)
class PacketCollection:
    def __init__(
        self,
        initial_radii,
        initial_nus,
        initial_mus,
        initial_energies,
        radiation_field_luminosity,
    ):
        self.initial_radii = initial_radii
        self.initial_nus = initial_nus
        self.initial_mus = initial_mus
        self.initial_energies = initial_energies

        self.radiation_field_luminosity = radiation_field_luminosity
        self.time_of_simulation = (
            1 / radiation_field_luminosity
        )  # 1 erg / luminosity
        self.output_nus = np.ones_like(initial_radii, dtype=np.float64) * -99.0
        self.output_energies = (
            np.ones_like(initial_radii, dtype=np.float64) * -99.0
        )


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
    ("last_interaction_shell_id", int64[:]),
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
        self.last_interaction_shell_id = -1 * np.ones(
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
        last_interaction_shell_id,
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
            temp_last_interaction_shell_id = np.empty(
                temp_length, dtype=np.int64
            )

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
            temp_last_interaction_shell_id[
                : self.length
            ] = self.last_interaction_shell_id

            self.nus = temp_nus
            self.energies = temp_energies
            self.initial_mus = temp_initial_mus
            self.initial_rs = temp_initial_rs
            self.last_interaction_in_nu = temp_last_interaction_in_nu
            self.last_interaction_type = temp_last_interaction_type
            self.last_interaction_in_id = temp_last_interaction_in_id
            self.last_interaction_out_id = temp_last_interaction_out_id
            self.last_interaction_shell_id = temp_last_interaction_shell_id
            self.length = temp_length

        self.nus[self.idx] = nu
        self.energies[self.idx] = energy
        self.initial_mus[self.idx] = initial_mu
        self.initial_rs[self.idx] = initial_r
        self.last_interaction_in_nu[self.idx] = last_interaction_in_nu
        self.last_interaction_type[self.idx] = last_interaction_type
        self.last_interaction_in_id[self.idx] = last_interaction_in_id
        self.last_interaction_out_id[self.idx] = last_interaction_out_id
        self.last_interaction_shell_id[self.idx] = last_interaction_shell_id
        self.idx += 1
