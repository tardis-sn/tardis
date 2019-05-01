from numba import jitclass, njit
import numpy as np
from astropy import constants as const

C_SPEED_OF_LIGHT = const.c.to('cm/s').value

@njit
def compute_distance2boundary(packet, storage):
    r = packet.r
    mu = packet.mu
    r_outer = storage.r_outer[packet.shell_id]
    r_inner = storage.r_inner[packet.shell_id]
    
    if (mu > 0.0):
        # direction outward 
        packet.shell_id += 1
        distance = np.sqrt(r_outer * r_outer + ((mu * mu - 1.0) * r * r)) - (r * mu)
    else:
        # going inward
        check = r_inner * r_inner + (r * r * (mu * my - 1.0))

        if (check >= 0.0):
            # hit inner boundary 
            packet.shell_id -= 1
            distance = -r * mu - np.sqrt(check)
        else:
            # miss inner boundary 
            packet.shell_id += 1
            distance = np.sqrt(r_outer * r_outer + ((mu * mu - 1.0) * r * r)) - (r * mu)
            
    packet.d_boundary = distance



# void
# compute_distance2boundary (rpacket_t * packet, const storage_model_t * storage)
# {
#   double r = rpacket_get_r (packet);
#   double mu = rpacket_get_mu (packet);
#   double r_outer = storage->r_outer[rpacket_get_current_shell_id (packet)];
#   double r_inner = storage->r_inner[rpacket_get_current_shell_id (packet)];
#   double check, distance;
#   if (mu > 0.0)
#     { // direction outward
#       rpacket_set_next_shell_id (packet, 1);
#       distance = sqrt (r_outer * r_outer + ((mu * mu - 1.0) * r * r)) - (r * mu);
#     }
#   else
#     { // going inward
#       if ( (check = r_inner * r_inner + (r * r * (mu * mu - 1.0)) )>= 0.0)
#         { // hit inner boundary
#           rpacket_set_next_shell_id (packet, -1);
#           distance = - r * mu - sqrt (check);
#         }
#       else
#         { // miss inner boundary
#           rpacket_set_next_shell_id (packet, 1);
#           distance = sqrt (r_outer * r_outer + ((mu * mu - 1.0) * r * r)) - (r * mu);
#         }
#     }
#   rpacket_set_d_boundary (packet, distance);
# }

# tardis_error_t


@njit
def compute_distance2line(rpacket, storage, c=C_SPEED_OF_LIGHT):
    if not rpacket_get_last_line(packet):
        r = rpacket.r
        mu = rpacket.mu
        nu = rpacket.nu
        nu_line = rpacket_get_nu_line(packet)
        
        distance = 0.0
        nu_diff = 0.0 

        ct =  storage.time_explosion * c
        doppler_factor = rpacket_doppler_factor(packet, storage)
        comov_nu = nu * doppler_factor

        nu_diff = comov_nu - nu_line
        if nu_diff >= 0:
            if not storage.full_relativity:
                distance = (nu_diff/nu) * ct
            else:
                nu_r = nu_line / nu
                distance = - mu * r + (ct - nu_r * nu_r * 
                    np.sqrt(ct * ct - (1 + r * r * (1 - mu * mu) * (1 + pow(nu_r, -2))))) / (1 + nu_r * nu_r)
            rpacket_set_d_line (packet, distance)
            return TARDIS_ERROR_OK
        else:
            if rpacket_get_next_line_id(packet) == storage.no_of_lines - 1:
                print("last_line = {}".format(storage.line_list_nu[rpacket_get_next_line_id(packet) - 1]))
                print("Last line in line list reached!")
            else if rpacket_get_next_line_id(packet) == 0:
                print("First line in line list!")
                print("next_line = {}".format(storage.line_list_nu[rpacket_get_next_line_id(packet) + 1]))
            else:
                print("last_line = {}".format(storage.line_list_nu[rpacket_get_next_line_id(packet) - 1]))
                print("next_line = {}".format(storage.line_list_nu[rpacket_get_next_line_id(packet) + 1]))
            print("ERROR: Comoving nu less than nu_line!")
            print("comov_nu = {}".format(comov_nu))
            print("nu_line = {}".format(nu_line))
            print("(comov_nu - nu_line) / nu_line = {}".format(comov_nu-nu_line/nu_line))
            print("r = {}".format(r))
            print("mu = {}".format(mu))
            print("nu = {}".format(nu))
            print("doppler_factor = {}".format(doppler_factor))
            print("cur_zone_id = {}".format(rpacket_get_current_shell_id(packet))
            return TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE
    else:
        rpacket_set_d_line(packet, MISS_DISTANCE)
        return TARDIS_ERROR_OK
            
                
            
# compute_distance2line (rpacket_t * packet, const storage_model_t * storage)
# {
#   if (!rpacket_get_last_line (packet))
#     {
#       double r = rpacket_get_r (packet);
#       double mu = rpacket_get_mu (packet);
#       double nu = rpacket_get_nu (packet);
#       double nu_line = rpacket_get_nu_line (packet);
#       double distance, nu_diff;
#       double ct = storage->time_explosion * c;
#       double doppler_factor = rpacket_doppler_factor (packet, storage);
#       double comov_nu = nu * doppler_factor;
#       if ( (nu_diff = comov_nu - nu_line) >= 0)
#         {
#           if (!storage->full_relativity)
#             {
#               distance = (nu_diff / nu) * ct;
#             }
#           else
#             {
#               double nu_r = nu_line / nu;
#               distance = - mu * r + (ct - nu_r * nu_r * sqrt(ct * ct -
#                 (1 + r * r * (1 - mu * mu) * (1 + pow(nu_r, -2))))) / (1 + nu_r * nu_r);
#             }
#           rpacket_set_d_line (packet, distance);
#           return TARDIS_ERROR_OK;
#         }
#       else
#         {
#           if (rpacket_get_next_line_id (packet) == storage->no_of_lines - 1)
#             {
#               fprintf (stderr, "last_line = %f\n",
#                        storage->
#                        line_list_nu[rpacket_get_next_line_id (packet) - 1]);
#               fprintf (stderr, "Last line in line list reached!");
#             }
#           else if (rpacket_get_next_line_id (packet) == 0)
#             {
#               fprintf (stderr, "First line in line list!");
#               fprintf (stderr, "next_line = %f\n",
#                        storage->
#                        line_list_nu[rpacket_get_next_line_id (packet) + 1]);
#             }
#           else
#             {
#               fprintf (stderr, "last_line = %f\n",
#                        storage->
#                        line_list_nu[rpacket_get_next_line_id (packet) - 1]);
#               fprintf (stderr, "next_line = %f\n",
#                        storage->
#                        line_list_nu[rpacket_get_next_line_id (packet) + 1]);
#             }
#           fprintf (stderr, "ERROR: Comoving nu less than nu_line!\n");
#           fprintf (stderr, "comov_nu = %f\n", comov_nu);
#           fprintf (stderr, "nu_line = %f\n", nu_line);
#           fprintf (stderr, "(comov_nu - nu_line) / nu_line = %f\n",
#                    (comov_nu - nu_line) / nu_line);
#           fprintf (stderr, "r = %f\n", r);
#           fprintf (stderr, "mu = %f\n", mu);
#           fprintf (stderr, "nu = %f\n", nu);
#           fprintf (stderr, "doppler_factor = %f\n", doppler_factor);
#           fprintf (stderr, "cur_zone_id = %" PRIi64 "\n", rpacket_get_current_shell_id (packet));
#           return TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE;
#         }
#     }
#   else
#     {
#       rpacket_set_d_line (packet, MISS_DISTANCE);
#       return TARDIS_ERROR_OK;
#     }
# }

@njit
def compute_distance2continuum(packet, storage):    
    # chi_continuum = 0.0
    # d_continuum = 0.0
    # chi_electron = storage.electron_densities[packet.shell_id] * storage.sigma_thomson

    # if (storage.full_relativity):
    #     chi_electron *= packet.doppler_factor(storage)
    
    # if (storage.cont_status == CONTINUUM_ON):
    #     if (packet.compute_chi_bf):
    #         calculate_chi_bf(packet, storage)
    #         calculate_chi_ff(packet, storage)
    #     else:
    #         packet.compute_chi_bf = True
    #     chi_continuum = packet.chi_boundfree + packet.chi_freefree + chi_electron 
    #     d_continuum = packet.tau_event / chi_continuum

    # else:
    #     chi_continuum = chi_electron
    #     d_continuum = storage.inverse_electron_densities[packet.shell_id] * \
    #         storage.inverse_sigma_thomson * packet.tau_event

    # if (packet.virtual_packet > 0):
    #     # set all continuum distances to MISS_DISTANCE in case of a virtual_packet 
    #     d_continuum = MISS_DISTANCE 
    #     packet.compute_chi_bf = False 
    # else:
    #     packet.chi_electron = chi_electron

    packet.d_continuum = storage.inverse_electron_densities[packet.shell_id] * \
            storage.inverse_sigma_thomson * packet.tau_event

# void
# compute_distance2continuum (rpacket_t * packet, storage_model_t * storage)
# {
#   double chi_continuum, d_continuum;
#   double chi_electron = storage->electron_densities[rpacket_get_current_shell_id(packet)] *
#     storage->sigma_thomson;
#   if (storage->full_relativity)
#     {
#       chi_electron *= rpacket_doppler_factor (packet, storage);
#     }

#   if (storage->cont_status == CONTINUUM_ON)
#     {
#     if (packet->compute_chi_bf)
#       {
#         calculate_chi_bf (packet, storage);
#         calculate_chi_ff (packet, storage);
#       }
#     else
#       {
#         packet->compute_chi_bf=true;
#       }
#       chi_continuum = rpacket_get_chi_boundfree (packet) + rpacket_get_chi_freefree (packet) + chi_electron;
#       d_continuum = rpacket_get_tau_event (packet) / chi_continuum;
#     }
#   else
#     {
#       chi_continuum = chi_electron;
#       d_continuum = storage->inverse_electron_densities[rpacket_get_current_shell_id (packet)] *
#         storage->inverse_sigma_thomson * rpacket_get_tau_event (packet);
#     }

#   if (rpacket_get_virtual_packet(packet) > 0)
#     {
#       //Set all continuum distances to MISS_DISTANCE in case of an virtual_packet
#       d_continuum = MISS_DISTANCE;
#       packet->compute_chi_bf = false;
#     }
#   else
#     {

#       //        fprintf(stderr, "--------\n");
#       //        fprintf(stderr, "nu = %e \n", rpacket_get_nu(packet));
#       //        fprintf(stderr, "chi_electron = %e\n", chi_electron);
#       //        fprintf(stderr, "chi_boundfree = %e\n", calculate_chi_bf(packet, storage));
#       //        fprintf(stderr, "chi_line = %e \n", rpacket_get_tau_event(packet) / rpacket_get_d_line(packet));
#       //        fprintf(stderr, "--------\n");

#       //rpacket_set_chi_freefree(packet, chi_freefree);
#       rpacket_set_chi_electron (packet, chi_electron);
#     }
#   rpacket_set_chi_continuum (packet, chi_continuum);
#   rpacket_set_d_continuum (packet, d_continuum);
# }