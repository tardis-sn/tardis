from numba import jitclass, njit
import numpy as np
from astropy import constants as const
from tardis.montecarlo.montecarlo_numba import njit_dict

C_SPEED_OF_LIGHT = const.c.to('cm/s').value
MISS_DISTANCE = 1e99

@njit(**njit_dict)
def compute_distance2boundary(packet, storage):
    r = packet.r
    mu = packet.mu
    r_outer = storage.r_outer[packet.current_shell_id]
    r_inner = storage.r_inner[packet.current_shell_id]
    
    if (mu > 0.0):
        # direction outward 
        packet.delta_shell_id = +1
        distance = np.sqrt(r_outer * r_outer + ((mu * mu - 1.0) * r * r)) - (r * mu)
    else:
        # going inward
        check = r_inner * r_inner + (r * r * (mu * mu - 1.0))

        if (check >= 0.0):
            # hit inner boundary 
            packet.delta_shell_id = -1
            distance = -r * mu - np.sqrt(check)
        else:
            # miss inner boundary 
            packet.delta_shell_id = + 1
            distance = np.sqrt(r_outer * r_outer + ((mu * mu - 1.0) * r * r)) - (r * mu)
            
    packet.d_boundary = distance



@njit
def compute_distance2line(packet, storage_model):
    if not packet.last_line:
        r = packet.r
        mu = packet.mu
        nu = packet.nu
        nu_line = packet.nu_line

        ct =  storage_model.time_explosion * C_SPEED_OF_LIGHT
        doppler_factor = packet.get_doppler_factor(storage_model)
        comov_nu = nu * doppler_factor

        nu_diff = comov_nu - nu_line
        if nu_diff >= 0:
            distance = (nu_diff/nu) * ct
            #else:
            #    nu_r = nu_line / nu
            #    distance = - mu * r + (ct - nu_r * nu_r * 
            #        np.sqrt(ct * ct - (1 + r * r * (1 - mu * mu) * (1 + pow(nu_r, -2))))) / (1 + nu_r * nu_r)
            packet.d_line = distance
        else:
            raise Exception
        #    if rpacket_get_next_line_id(packet) == storage.no_of_lines - 1:
        #        print("last_line = {}".format(storage.line_list_nu[rpacket_get_next_line_id(packet) - 1]))
        #        print("Last line in line list reached!")
        #    elif rpacket_get_next_line_id(packet) == 0:
        #        print("First line in line list!")
        #        print("next_line = {}".format(storage.line_list_nu[rpacket_get_next_line_id(packet) + 1]))
        #    else:
        #        print("last_line = {}".format(storage.line_list_nu[rpacket_get_next_line_id(packet) - 1]))
        #        print("next_line = {}".format(storage.line_list_nu[rpacket_get_next_line_id(packet) + 1]))
        #    print("ERROR: Comoving nu less than nu_line!")
        #    print("comov_nu = {}".format(comov_nu))
        #    print("nu_line = {}".format(nu_line))
        #    print("(comov_nu - nu_line) / nu_line = {}".format(comov_nu-nu_line/nu_line))
        #    print("r = {}".format(r))
        #    print("mu = {}".format(mu))
        #    print("nu = {}".format(nu))
        #    print("doppler_factor = {}".format(doppler_factor))
        #    print("cur_zone_id = {}".format(rpacket_get_current_shell_id(packet))
        #    #return TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE
    else:
        packet.d_line = MISS_DISTANCE
        #return TARDIS_ERROR_OK

@njit(**njit_dict)
def compute_distance2continuum(packet, storage):    
    packet.d_electron = storage.inverse_electron_densities[packet.current_shell_id] * \
            storage.inverse_sigma_thomson * packet.tau_event

