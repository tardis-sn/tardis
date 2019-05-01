import numpy as np
from enum import Enum
from numba import int64, float64
from numba import jitclass, njit
from tardis.montecarlo.montecarlo_numba.compute_distance import (compute_distance2boundary,
         
                                                                 compute_distance2continuum)

from astropy import constants as const



C_SPEED_OF_LIGHT = const.c.to('cm/s').value

class PacketStatus(Enum):
  IN_PROCESS = 0
  EMITTED = 1
  REABSORBED = 2

class InteractionType(Enum):
    ESCATTERING = 0
    BOUNDARY = 1
    LINE = 2


rpacket_spec = [
    ('r', float64),
    ('mu', float64),
    ('nu', float64),
    ('energy', float64)
    ('tau_event', float64),
 #   ('nu_line', float64)
 #   ('d_line', float64)  # Distance to electron event. 
    ('d_electron', float64), #/**< Distance to line event. */
    ('d_boundary', float64), # distance to boundary 
    ('current_shell_id', int64),
    ('next_shell_id', int64),
    ('next_interaction', enum)
]


@jitclass(rpacket_spec)
class RPacket(object):
    def __init__(self, r, mu, nu):
        self.r = r
        self.mu = mu
        self.nu = nu
        self.energy = energy
        self.tau_event = self.get_tau_event()
        #self.nu_line = line_search() # TODO: Implement this
        self.current_shell_id = 0
        self.compute_distances()
        #self.nu_line = nu_line

    @staticmethod
    def get_tau_event():
        return  np.random.exponential()
    
    def compute_distances(self, storage):
        compute_distance2continuum(self, storage)
        compute_distance2boundary(self, storage)
        if packet.d_boundary < packet.d_electron:
            next_interaction = InteractionType.BOUNDARY
        else:
            next_interaction = InteractionType.ESCATTERING
        packet.next_interaction = next_interaction

#     double
# rpacket_doppler_factor (const rpacket_t *packet, const storage_model_t *storage)
# {
#   double beta = rpacket_get_r (packet) * storage->inverse_time_explosion * INVERSE_C;
#   if (!storage->full_relativity)
#     {
#       return 1.0 - rpacket_get_mu (packet) * beta;
#     }
#   else
#     {
#       return (1.0 - rpacket_get_mu (packet) * beta) / sqrt (1 - beta * beta);
#     }
# }
    
    def get_doppler_factor(self, storage):
        beta = self.r * storage.inverse_time_explosion / C_SPEED_OF_LIGHT
        
        if not storage.full_relativity:
            return 1.0 - self.mu * beta
        else:
            return (1.0 - self.mu * beta) / np.sqrt(1 - beta*beta)

    def move_packet(self, storage, distance):
        doppler_factor = self.get_doppler_factor(storage)
        r = packet.r
        if (distance > 0.0):
            new_r = np.sqrt (r * r + distance * distance +
                             2.0 * r * distance * packet.mu)
            packet.mu = (packet.mu * r + distance) / new_r
            packet.r = new_r


    def move_packet_across_shell_boundary(self, storage, distance, mt_state):
        self.move_packet(storage, distance)

        if self.virtual_packet > 0:
            delta_tau_event = self.chi_continuum * distance
            self.tau_event = self.tau_event + delta_tau_event
            self.compute_chi_bf = True
        else:
            self.reset_tau_event(mt_state)

        if ((self.current_shell_id < storage.no_of_shells - 1) and packet.next_shell_id == 1) or ()

    void
move_packet_across_shell_boundary (rpacket_t * packet,
                                   storage_model_t * storage, double distance, rk_state *mt_state)
{
  move_packet (packet, storage, distance);
  if (rpacket_get_virtual_packet (packet) > 0)
    {
      double delta_tau_event = rpacket_get_chi_continuum(packet) * distance;
      rpacket_set_tau_event (packet,
                             rpacket_get_tau_event (packet) +
                             delta_tau_event);
	  packet->compute_chi_bf = true;
    }
  else
    {
      rpacket_reset_tau_event (packet, mt_state);
    }
  if ((rpacket_get_current_shell_id (packet) < storage->no_of_shells - 1
       && rpacket_get_next_shell_id (packet) == 1)
      || (rpacket_get_current_shell_id (packet) > 0
          && rpacket_get_next_shell_id (packet) == -1))
    {
      rpacket_set_current_shell_id (packet,
                                    rpacket_get_current_shell_id (packet) +
                                    rpacket_get_next_shell_id (packet));
    }
  else if (rpacket_get_next_shell_id (packet) == 1)
    {
      rpacket_set_status (packet, TARDIS_PACKET_STATUS_EMITTED);
    }
  else if ((storage->reflective_inner_boundary == 0) ||
           (rk_double (mt_state) > storage->inner_boundary_albedo))
    {
      rpacket_set_status (packet, TARDIS_PACKET_STATUS_REABSORBED);
    }
  else
    {
      double doppler_factor = rpacket_doppler_factor (packet, storage);
      double comov_nu = rpacket_get_nu (packet) * doppler_factor;
      double comov_energy = rpacket_get_energy (packet) * doppler_factor;
      // TODO: correct
      rpacket_set_mu (packet, rk_double (mt_state));
      double inverse_doppler_factor = rpacket_inverse_doppler_factor (packet, storage);
      rpacket_set_nu (packet, comov_nu * inverse_doppler_factor);
      rpacket_set_energy (packet, comov_energy * inverse_doppler_factor);
      if (rpacket_get_virtual_packet_flag (packet) > 0)
        {
          montecarlo_one_packet (storage, packet, -2, mt_state);
        }
    }
}