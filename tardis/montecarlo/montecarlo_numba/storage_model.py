from numba import int64, float64, jitclass
import numpy as np
from astropy import constants as const

C_SPEED_OF_LIGHT = const.c.to('cm/s').value

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
    ('electron_densities', float64[:]),
    ('inverse_electron_densities', float64[:]), # Maybe remove the inverse things
    ('line_list_nu', float64[:]),
    ('line_lists_tau_sobolevs', float64[:, :]),
    ('no_of_lines', int64),
    ('line_interaction_id', int64),
#    ('*js', float64),
#    ('*nubars', float64),
    ('sigma_thomson', float64),
    ('inverse_sigma_thomson', float64),
    ('ct', float64),
]

@jitclass(storage_model_spec)
class StorageModel(object):
    def __init__(self, packet_nus, packet_mus, packet_energies, 
    output_nus, output_energies, no_of_packets, no_of_shells, 
    r_inner, r_outer, v_inner, time_explosion, electron_densities, 
    line_list_nu, line_lists_tau_sobolevs, no_of_lines, line_interaction_id, sigma_thomson):
        self.packet_nus = packet_nus
        self.packet_mus = packet_mus
        self.packet_energies = packet_energies
        self.no_of_packets = len(self.packet_nus)
        
        self.output_nus = output_nus
        self.output_energies = output_energies
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.v_inner = v_inner
        self.no_of_shells = len(self.v_inner)

        self.time_explosion = time_explosion
        self.inverse_time_explosion = 1 / time_explosion

        self.electron_densities = electron_densities

        self.inverse_electron_densities = 1 / electron_densities
        
        self.sigma_thomson = sigma_thomson

        self.inverse_sigma_thomson = 1 / self.sigma_thomson
        self.no_of_lines = no_of_lines
        self.line_list_nu = line_list_nu
        self.line_lists_tau_sobolevs = line_lists_tau_sobolevs
        self.ct = self.time_explosion * C_SPEED_OF_LIGHT

def initialize_storage_model(model, plasma, runner):
    storage_model_kwargs = {'packet_nus': runner.input_nu,
    'packet_mus': runner.input_mu,
    'packet_energies': runner.input_energy,
    'output_nus': runner._output_nu,
    'output_energies': runner._output_energy,
    'no_of_packets': runner.input_nu.size,
    'no_of_shells': model.no_of_shells,
    'r_inner': runner.r_inner_cgs,
    'r_outer': runner.r_outer_cgs,
    'v_inner': runner.v_inner_cgs,
    'time_explosion': model.time_explosion.to('s').value,
    'electron_densities': plasma.electron_densities.values,
    'line_list_nu': plasma.atomic_data.lines.nu.values, 
    'no_of_lines': len(plasma.atomic_data.lines.nu.values),
    'line_interaction_id': runner.get_line_interaction_id(
        runner.line_interaction_type),
    'sigma_thomson': runner.sigma_thomson.cgs.value}
    storage_model_kwargs['line_lists_tau_sobolevs']= np.ascontiguousarray(
        plasma.tau_sobolevs.values.copy(), dtype=np.float64)
    return StorageModel(**storage_model_kwargs)