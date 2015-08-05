from astropy import units as u
from astropy.utils import lazyproperty

from tardis.montecarlo import montecarlo

import numpy as np

class MontecarloRunner(object):
    """
    This class is designed as an interface between the Python part and the
    montecarlo C-part
    """

    def run(self, model, no_of_virtual_packets, nthreads=1):
        self.time_of_simulation = model.time_of_simulation
        montecarlo.montecarlo_radial1d(
            model, self, virtual_packet_flag=no_of_virtual_packets,
            nthreads=nthreads)

    def legacy_return(self):
        return (self.packet_nu, self.packet_energy,
                self.j_estimator, self.nu_bar_estimator,
                self.last_line_interaction_in_id,
                self.last_line_interaction_out_id,
                self.last_interaction_type,
                self.last_line_interaction_shell_id)

    @property
    def packet_nu(self):
        return u.Quantity(self._packet_nu, u.Hz)

    @property
    def packet_energy(self):
        return u.Quantity(self._packet_energy, u.erg)

    @property
    def packet_luminosity(self):
        return self.packet_energy / self.time_of_simulation

    @property
    def emitted_packet_mask(self):
        return self.packet_energy >=0

    @property
    def emitted_packet_nu(self):
        return self.packet_nu[self.emitted_packet_mask]

    @property
    def reabsorbed_packet_nu(self):
        return self.packet_nu[~self.emitted_packet_mask]

    @property
    def reabsorbed_packet_luminosity(self):
        return -self.packet_luminosity[~self.emitted_packet_mask]


    @property
    def emitted_packet_luminosity(self):
        return self.packet_luminosity[self.emitted_packet_mask]

    @property
    def reabsorbed_packet_luminosity(self):
        return -self.packet_luminosity[~self.emitted_packet_mask]


    @staticmethod
    def generate_spectrum(nu, energy, nu_bins):
        return np.histogram(nu, weights=energy, bins=nu_bins)