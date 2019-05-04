vpacket_spec = [
    ('r', float64),
    ('mu', float64),
    ('nu', float64),
    ('energy', float64),
    ('next_line_id', int64),
    ('current_shell_id', int64),
    ('status', int64),
]

@jitclass(vpacket_spec)
class VPacket(object):
    def __init__(self, r, mu, nu, energy):
        self.r = r
        self.mu = mu
        self.nu = nu
        self.energy = energy
        self.current_shell_id = 0
        self.status = IN_PROCESS
        self.next_line_id = -1

    def trace_packet(self, storage_model):
        
        r_inner = storage_model.r_inner[self.current_shell_id]
        r_outer = storage_model.r_outer[self.current_shell_id]
        
        distance = 0.0

        distance_boundary, delta_shell = calculate_distance_boundary(
            self.r, self.mu, r_inner, r_outer)
        
        #defining start for line interaction
        cur_line_id = self.next_line_id
        nu_line = 0.0

        #defining taus
        tau_event = np.random.exponential()
        tau_trace_line = 0.0
        tau_trace_line_combined = 0.0
        
        #e scattering initialization

        cur_electron_density = storage_model.electron_densities[
            self.current_shell_id]
        cur_inverse_electron_density = 1 / cur_electron_density
        distance_electron = calculate_distance_electron(
            cur_inverse_electron_density, tau_event)


        #Calculating doppler factor
        doppler_factor = get_doppler_factor(self.r, self.mu, 
                                        storage_model.inverse_time_explosion)
        comov_nu = self.nu * doppler_factor
        distance_trace = 0.0
        last_line = False

        while True:
            if cur_line_id < storage_model.no_of_lines: # not last_line
                nu_line = storage_model.line_list_nu[cur_line_id]
                tau_trace_line = storage_model.line_lists_tau_sobolevs[cur_line_id, 
                        self.current_shell_id]
            else:
                last_line = True
                self.next_line_id = cur_line_id
                break
            
            tau_trace_line_combined += tau_trace_line
            distance_trace = calculate_distance_line(self.nu, comov_nu, nu_line, 
                                                        storage_model.ct)
            tau_trace_electron = calculate_tau_electron(cur_electron_density, 
                                                        distance_trace)

            tau_trace_combined = tau_trace_line_combined + tau_trace_electron

            if (distance_boundary <= distance_trace):
                ## current_shell_id +=1
                ## distance_boundary
                #unless shell 
                interaction_type = InteractionType.BOUNDARY # BOUNDARY
                self.next_line_id = cur_line_id
                distance = distance_boundary
                break
            
            if (distance_electron < distance_trace) and (distance_electron < distance_boundary):
                interaction_type = InteractionType.ESCATTERING
                distance = distance_electron
                self.next_line_id = cur_line_id
                break
                        
            cur_line_id += 1
        if not last_line:            
            return distance, interaction_type, delta_shell
        else:
            if distance_electron < distance_boundary:
                return distance_electron, InteractionType.ESCATTERING, delta_shell
            else:
                return distance_boundary, InteractionType.BOUNDARY, delta_shell
