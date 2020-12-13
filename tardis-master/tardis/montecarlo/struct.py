from ctypes import Structure, POINTER, c_int, c_int64, c_double, c_ulong, c_bool
from enum import Enum

c_tardis_error_t = c_int
c_rpacket_status_t = c_int
c_cont_status_t = c_int


class RPacket(Structure):
    _fields_ = [
        ("nu", c_double),
        ("mu", c_double),
        ("energy", c_double),
        ("r", c_double),
        ("tau_event", c_double),
        ("nu_line", c_double),
        ("current_shell_id", c_int64),
        ("next_line_id", c_int64),
        ("last_line", c_int64),
        ("close_line", c_int64),
        ("current_continuum_id", c_int64),
        ("virtual_packet_flag", c_int64),
        ("virtual_packet", c_int64),
        ("d_line", c_double),
        ("d_electron", c_double),
        ("d_boundary", c_double),
        ("d_cont", c_double),
        ("next_shell_id", c_int64),
        ("status", c_rpacket_status_t),
        ("id", c_int64),
        ("chi_th", c_double),
        ("chi_cont", c_double),
        ("chi_ff", c_double),
        ("chi_bf", c_double),
        ("chi_bf_tmp_partial", POINTER(c_double)),
        ("macro_atom_activation_level", c_int64),
        ("compute_chi_bf", c_bool),
        ("vpacket_weight", c_double),
    ]


class PhotoXsect1level(Structure):
    _fields_ = [
        ("nu", POINTER(c_double)),
        ("x_sect", POINTER(c_double)),
        ("no_of_points", c_int64),
    ]


class StorageModel(Structure):
    _fields_ = [
        ("packet_nus", POINTER(c_double)),
        ("packet_mus", POINTER(c_double)),
        ("packet_energies", POINTER(c_double)),
        ("output_nus", POINTER(c_double)),
        ("output_energies", POINTER(c_double)),
        ("last_interaction_in_nu", POINTER(c_double)),
        ("last_line_interaction_in_id", POINTER(c_int64)),
        ("last_line_interaction_out_id", POINTER(c_int64)),
        ("last_line_interaction_shell_id", POINTER(c_int64)),
        ("last_interaction_type", POINTER(c_int64)),
        ("last_interaction_out_type", POINTER(c_int64)),
        ("no_of_packets", c_int64),
        ("no_of_shells", c_int64),
        ("no_of_shells_i", c_int64),
        ("r_inner", POINTER(c_double)),
        ("r_outer", POINTER(c_double)),
        ("r_inner_i", POINTER(c_double)),
        ("r_outer_i", POINTER(c_double)),
        ("v_inner", POINTER(c_double)),
        ("time_explosion", c_double),
        ("inverse_time_explosion", c_double),
        ("electron_densities", POINTER(c_double)),
        ("electron_densities_i", POINTER(c_double)),
        ("inverse_electron_densities", POINTER(c_double)),
        ("line_list_nu", POINTER(c_double)),
        ("continuum_list_nu", POINTER(c_double)),
        ("line_lists_tau_sobolevs", POINTER(c_double)),
        ("line_lists_tau_sobolevs_i", POINTER(c_double)),
        ("line_lists_tau_sobolevs_nd", c_int64),
        ("line_lists_j_blues", POINTER(c_double)),
        ("line_lists_j_blues_nd", c_int64),
        ("line_lists_Edotlu", POINTER(c_double)),
        ("no_of_lines", c_int64),
        ("no_of_edges", c_int64),
        ("line_interaction_id", c_int64),
        ("transition_probabilities", POINTER(c_double)),
        ("transition_probabilities_nd", c_int64),
        ("line2macro_level_upper", POINTER(c_int64)),
        ("macro_block_references", POINTER(c_int64)),
        ("transition_type", POINTER(c_int64)),
        ("destination_level_id", POINTER(c_int64)),
        ("transition_line_id", POINTER(c_int64)),
        ("js", POINTER(c_double)),
        ("nubars", POINTER(c_double)),
        ("spectrum_start_nu", c_double),
        ("spectrum_delta_nu", c_double),
        ("spectrum_end_nu", c_double),
        ("spectrum_virt_start_nu", c_double),
        ("spectrum_virt_end_nu", c_double),
        ("spectrum_virt_nu", POINTER(c_double)),
        ("sigma_thomson", c_double),
        ("inverse_sigma_thomson", c_double),
        ("inner_boundary_albedo", c_double),
        ("reflective_inner_boundary", c_int64),
        ("current_packet_id", c_int64),
        ("photo_xsect", POINTER(POINTER(PhotoXsect1level))),
        ("chi_ff_factor", POINTER(c_double)),
        ("t_electrons", POINTER(c_double)),
        ("l_pop", POINTER(c_double)),
        ("l_pop_r", POINTER(c_double)),
        ("cont_status", c_cont_status_t),
        ("bf_treatment", c_int),
        ("virt_packet_nus", POINTER(c_double)),
        ("virt_packet_energies", POINTER(c_double)),
        ("virt_packet_last_interaction_in_nu", POINTER(c_double)),
        ("virt_packet_last_interaction_type", POINTER(c_int64)),
        ("virt_packet_last_line_interaction_in_id", POINTER(c_int64)),
        ("virt_packet_last_line_interaction_out_id", POINTER(c_int64)),
        ("virt_packet_count", c_int64),
        ("virt_array_size", c_int64),
        ("kpacket2macro_level", c_int64),
        ("cont_edge2macro_level", POINTER(c_int64)),
        ("photo_ion_estimator", POINTER(c_double)),
        ("stim_recomb_estimator", POINTER(c_double)),
        ("photo_ion_estimator_statistics", POINTER(c_int64)),
        ("bf_heating_estimator", POINTER(c_double)),
        ("ff_heating_estimator", POINTER(c_double)),
        ("stim_recomb_cooling_estimator", POINTER(c_double)),
        ("full_relativity", c_int),
        ("survival_probability", c_double),
        ("tau_russian", c_double),
        ("tau_bias", POINTER(c_double)),
        ("enable_biasing", c_int),
    ]


# Variables corresponding to `tardis_error_t` enum.
TARDIS_ERROR_OK = 0
TARDIS_ERROR_BOUNDS_ERROR = 1
TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE = 2

# Variables corresponding to `rpacket_status_t` enum.
TARDIS_PACKET_STATUS_IN_PROCESS = 0
TARDIS_PACKET_STATUS_EMITTED = 1
TARDIS_PACKET_STATUS_REABSORBED = 2

# Variables corresponding to `ContinuumProcessesStatus` enum.
CONTINUUM_OFF = 0
CONTINUUM_ON = 1

# Variables corresponding to `emission_type` enum.
BB_EMISSION = -1
BF_EMISSION = -2
FF_EMISSION = -3

# Variables corresponding to `bound_free_treatment` enum.
class BoundFreeTreatment(Enum):
    LIN_INTERPOLATION = 0
    HYDROGENIC = 1


# Variables corresponding to macros defined in rpacket.h .
MISS_DISTANCE = 1e99
C = 29979245800.0
INVERSE_C = 3.33564095198152e-11
H = 6.6260755e-27
KB = 1.3806488e-16

# `rk_state` specific macros.
RK_STATE_LEN = 624


class RKState(Structure):
    _fields_ = [
        ("key", c_ulong * RK_STATE_LEN),
        ("pos", c_int),
        ("has_gauss", c_int),
        ("gauss", c_double),
    ]
