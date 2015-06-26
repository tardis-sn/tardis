from tardis.montecarlo import montecarlo

class MontecarloRunner():

    def __init__(self, no_of_virtual_packets, tardis_config, packet_src,
        plasma_array, atom_data, j_blue_estimators,
        montecarlo_virtual_luminosity, transition_probabilities=None):
        self.no_of_virtual_packets = no_of_virtual_packets
        self.tardis_config = tardis_config
        self.packet_src = packet_src
        self.plasma_array = plasma_array
        self.atom_data = atom_data
        self.j_blue_estimators = j_blue_estimators
        self.montecarlo_virtual_luminosity = montecarlo_virtual_luminosity
        self.plasma_array.tau_sobolevs.update(self.plasma_array.tau_sobolevs.copy('F'))
        self.plasma_array.electron_densities.update(self.plasma_array.electron_densities.copy('F'))
        if transition_probabilities is not None:
            self.transition_probabilities = transition_probabilities.copy('F')

    def run_montecarlo_radial1d(self):
        montecarlo_nu, montecarlo_energies, j_estimators, nubar_estimators, \
        last_line_interaction_in_id, last_line_interaction_out_id, \
        last_interaction_type, last_line_interaction_shell_id = \
        montecarlo.montecarlo_radial1d(self, virtual_packet_flag=
            self.no_of_virtual_packets)
        return montecarlo_nu, montecarlo_energies, j_estimators, \
            nubar_estimators, last_line_interaction_in_id,\
            last_line_interaction_out_id, last_interaction_type,\
            last_line_interaction_shell_id


