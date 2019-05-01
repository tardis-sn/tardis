from numba import int64, float64
from astropy import constants as const
storage_model_spec = [
    ('packet_nus', float64[:]),
    ('packet_mus', float64[:]),
    ('packet_energies', float64[:]),
    ('output_nus', float64[:]),
    ('output_energies', float64[:]),
    ('no_of_packets', int64),
    ('no_of_shells', int64),
    ('r_inner', float64[:]),
    ('r_outer', float64[:]),
    ('v_inner', float64[:]),
    ('time_explosion', float64),
    ('inverse_time_explosion', float64),
    ('*electron_densities', float64),
    ('*inverse_electron_densities', float64),
    ('*line_list_nu', float64),
    ('*line_lists_tau_sobolevs', float64),
    ('line_lists_tau_sobolevs_nd', int64),
    ('no_of_lines', int64),
    ('no_of_edges', int64),
    ('line_interaction_id', int64),
#    ('*js', float64),
#    ('*nubars', float64),
    ('sigma_thomson', float64),
    ('inverse_sigma_thomson', float64),
]

class StorageModel(object):
    def __init__(self, packet_nus, packet_mus, packet_energies, 
    output_nus, output_energies, no_of_packets, no_of_shells, 
    r_inner, r_outer, v_inner, time_explosion, electron_densities, line_list_nu, line_lists_tau_sobolevs, line_lists_tau_sobolevs_nd, 
    no_of_lines, no_of_edges, line_interaction_id, 
    inverse_sigma_thomson):
        self.packet_nus = packet_nus
        self.packet_mus = packet_mus
        self.packet_energies = packet_energies
        self.output_nus = output_nus
        self.output_energies = output_energies
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.v_inner = v_inner
        
        self.time_explosion = time_explosion
        self.inverse_time_explosion = 1 / time_explosion

        self.electron_densities = electron_densities

        self.inverse_electron_densities = 1 / electron_densities
        
        self.sigma_thomson = const.sigma_T.to('cm^2').value
        self.inverse_sigma_thomson = 1 / self.sigma_thomson
