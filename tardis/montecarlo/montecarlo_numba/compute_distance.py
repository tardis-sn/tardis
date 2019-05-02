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
    if distance < 0.0:
        print('hello')
        pass
    packet.d_boundary = distance



@njit
def compute_distance2line(packet, storage_model):
    if not packet.last_line and not packet.close_line:
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
    else:
        packet.d_line = MISS_DISTANCE
        #return TARDIS_ERROR_OK

@njit(**njit_dict)
def compute_distance2continuum(packet, storage):    
    packet.d_electron = storage.inverse_electron_densities[packet.current_shell_id] * \
            storage.inverse_sigma_thomson * packet.tau_event

