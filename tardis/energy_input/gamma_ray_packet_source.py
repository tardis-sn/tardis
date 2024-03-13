import numpy as np
from numba import njit

from tardis.energy_input.energy_source import (
    positronium_continuum,
)
from tardis.energy_input.GXPacket import (
    GXPacketCollection,
    initialize_packet_properties,
)
from tardis.energy_input.util import (
    get_index,
    get_random_unit_vector,
    doppler_factor_3d,
    H_CGS_KEV,
)
from tardis.energy_input.samplers import sample_energy
from tardis.montecarlo.montecarlo_numba import njit_dict_no_parallel
from tardis.montecarlo.packet_source import BasePacketSource


class RadioactivePacketSource(BasePacketSource):
    def __init__(
        self,
        packet_energy,
        gamma_ray_lines,
        positronium_fraction,
        inner_velocities,
        outer_velocities,
        inv_volume_time,
        times,
        energy_df_rows,
        effective_times,
        taus,
        parents,
        average_positron_energies,
        inventories,
        average_power_per_mass,
        **kwargs,
    ):
        self.packet_energy = packet_energy
        self.gamma_ray_lines = gamma_ray_lines
        self.positronium_fraction = positronium_fraction
        self.inner_velocities = inner_velocities
        self.outer_velocities = outer_velocities
        self.inv_volume_time = inv_volume_time
        self.times = times
        self.energy_df_rows = energy_df_rows
        self.effective_times = effective_times
        self.taus = taus
        self.parents = parents
        self.average_positron_energies = average_positron_energies
        self.inventories = inventories
        self.average_power_per_mass = average_power_per_mass
        super().__init__(**kwargs)

    @njit(**njit_dict_no_parallel)
    def create_packet_radii(
        self, no_of_packets, inner_velocity, outer_velocity
    ):
        """Initialize the random radii of packets in a shell

        Parameters
        ----------
        packet_count : int
            Number of packets in the shell
        inner_velocity : float
            Inner velocity of the shell
        outer_velocity : float
            Outer velocity of the shell

        Returns
        -------
        array
            Array of length packet_count of random locations in the shell
        """
        z = np.random.random(no_of_packets)
        initial_radii = (
            z * inner_velocity**3.0 + (1.0 - z) * outer_velocity**3.0
        ) ** (1.0 / 3.0)

        return initial_radii

    def create_packet_nus(
        self,
        no_of_packets,
        energy,
        intensity,
        positronium_fraction,
        positronium_energy,
        positronium_intensity,
    ):
        """Create an array of packet frequency-energies (i.e. E = h * nu)

        Parameters
        ----------
        no_of_packets : int
            Number of packets to produce frequency-energies for
        energy : One-dimensional Numpy Array, dtype float
            Array of frequency-energies to sample
        intensity : One-dimensional Numpy Array, dtype float
            Array of intensities to sample
        positronium_fraction : float
            The fraction of positrons that form positronium
        positronium_energy : array
            Array of positronium frequency-energies to sample
        positronium_intensity : array
            Array of positronium intensities to sample

        Returns
        -------
        array
            Array of sampled frequency-energies
        array
            Positron creation mask
        """
        nu_energies = np.zeros(no_of_packets)
        positrons = np.zeroes(no_of_packets)
        zs = np.random.random(no_of_packets)
        for i in range(no_of_packets):
            nu_energies[i] = sample_energy(energy, intensity)
            # positron
            if nu_energies[i] == 511:
                # positronium formation 25% of the time if fraction is 1
                if zs[i] < positronium_fraction and np.random.random() < 0.25:
                    nu_energies[i] = sample_energy(
                        positronium_energy, positronium_intensity
                    )
                positrons[i] = 1

        return nu_energies, positrons

    def create_packet_directions(self, no_of_packets):
        """Create an array of random directions

        Parameters
        ----------
        no_of_packets : int
            Number of packets to produce directions for

        Returns
        -------
        array
            Array of direction vectors
        """
        directions = np.zeros((no_of_packets, 3))
        for i in range(no_of_packets):
            directions[i, :] = get_random_unit_vector()

        return directions

    def create_packet_energies(self, no_of_packets, energy):
        """Create the uniform packet energy for a number of packets

        Parameters
        ----------
        no_of_packets : int
            Number of packets
        energy : float
            The packet energy

        Returns
        -------
        array
            Array of packet energies
        """
        return np.ones(no_of_packets) * energy

    def create_packet_times_uniform_time(self, no_of_packets, start, end):
        """Samples decay time uniformly (needs non-uniform packet energies)

        Parameters
        ----------
        no_of_packets : int
            Number of packets
        start : float
            Start time
        end : float
            End time

        Returns
        -------
        array
            Array of packet decay times
        """
        z = np.random.random(no_of_packets)
        decay_times = z * start + (1 - z) * end
        return decay_times

    @njit(**njit_dict_no_parallel)
    def create_packet_times_uniform_energy(
        self,
        no_of_packets,
        start_tau,
        end_tau=0.0,
        decay_time_min=0.0,
        decay_time_max=0.0,
    ):
        """Samples the decay time from the mean lifetime of the isotopes

        Parameters
        ----------
        no_of_packets : int
            Number of packets
        start_tau : float
            Initial isotope mean lifetime
        end_tau : float, optional
            Ending mean lifetime, by default 0.0 for single decays
        decay_time_min : float, optional
            Minimum time to decay, by default 0.0
        decay_time_max : float, optional
            Maximum time to decay, by default 0.0

        Returns
        -------
        array
            Array of decay times
        """
        decay_times = np.ones(no_of_packets) * decay_time_min
        for i in range(no_of_packets):
            # rejection sampling
            while (decay_times[i] <= decay_time_min) or (
                decay_times[i] >= decay_time_max
            ):
                decay_times[i] = -start_tau * np.log(
                    np.random.random()
                ) - end_tau * np.log(np.random.random())
        return decay_times

    @njit(**njit_dict_no_parallel)
    def calculate_energy_factors(self, no_of_packets, start_time, decay_times):
        energy_factors = np.ones(no_of_packets)
        for i in range(no_of_packets):
            if decay_times[i] < start_time:
                energy_factors[i] = decay_times[i] / start_time
                decay_times[i] = start_time
        return energy_factors, decay_times

    def create_packets(self, decays_per_isotope, *args, **kwargs):
        """Initialize a collection of GXPacket objects for the simulation
        to operate on.

        Parameters
        ----------
        decays_per_isotope : array int64
            Number of decays per simulation shell per isotope

        Returns
        -------
        list
            List of GXPacket objects
        array
            Array of main output dataframe rows
        array
            Array of plotting output dataframe rows
        array
            Array of positron output dataframe rows
        """
        number_of_packets = decays_per_isotope.sum().sum()
        decays_per_shell = decays_per_isotope.T.sum().values

        energy_plot_df_rows = np.zeros((number_of_packets, 8))
        energy_plot_positron_rows = np.zeros((number_of_packets, 4))

        positronium_energy, positronium_intensity = positronium_continuum()

        packet_index = 0
        # go through each shell
        for k, shell in enumerate(decays_per_shell):
            isotope_packet_count_df = decays_per_isotope.iloc[k]

            # packet index
            i = 0
            for (
                isotope_name,
                isotope_packet_count,
            ) in isotope_packet_count_df.items():
                # get isotope information
                isotope_energy = self.gamma_ray_lines[isotope_name][0, :]
                isotope_intensity = self.gamma_ray_lines[isotope_name][1, :]
                isotope_positron_fraction = self.calculate_positron_fraction(
                    self.average_positron_energies[isotope_name],
                    isotope_energy,
                    isotope_intensity,
                )
                tau_start = self.taus[isotope_name]

                if isotope_name in self.parents:
                    tau_end = self.taus[self.parents[isotope_name]]
                else:
                    tau_end = 0

                # sample radii at time = 0
                initial_radii = self.create_packet_radii(
                    isotope_packet_count,
                    self.inner_velocities[k],
                    self.outer_velocities[k],
                )

                # sample directions (valid at all times)
                initial_directions = self.create_packet_directions(
                    isotope_packet_count
                )

                # packet decay time
                initial_times = self.create_packet_times_uniform_energy(
                    isotope_packet_count,
                    tau_start,
                    tau_end,
                    decay_time_min=0,
                    decay_time_max=self.times[-1],
                )

                # get the time step index of the packets
                initial_time_indexes = np.array(
                    [
                        get_index(decay_time, self.times)
                        for decay_time in initial_times
                    ]
                )

                # get the time of the middle of the step for each packet
                packet_effective_times = np.array(
                    [self.effective_times[i] for i in initial_time_indexes]
                )

                # scale radius by packet decay time. This could be replaced with
                # Geometry object calculations. Note that this also adds a random
                # unit vector multiplication for 3D. May not be needed.
                initial_locations = (
                    initial_radii
                    * packet_effective_times
                    * self.create_packet_directions(isotope_packet_count)
                )

                # get the packet shell index
                initial_shells = np.ones(isotope_packet_count) * k

                # the individual gamma-ray energies that make up a packet
                # co-moving frame, including positronium formation
                initial_nu_energies_cmf, positron_mask = self.create_packet_nus(
                    isotope_packet_count,
                    isotope_energy,
                    isotope_intensity,
                    self.positronium_fraction,
                    positronium_energy,
                    positronium_intensity,
                )

                # equivalent frequencies
                initial_nus_cmf = initial_nu_energies_cmf / H_CGS_KEV

                # compute scaling factor for packets emitted before start time
                # and move packets to start at that time
                # probably not necessary- we have rejection sampling in the
                # create_packet_times_uniform_energy method
                energy_factors, initial_times = self.calculate_energy_factors(
                    isotope_packet_count, self.times[0], initial_times
                )

                # the CMF energy of a packet scaled by the "early energy factor"
                initial_packet_energies_cmf = (
                    self.create_packet_energies(
                        isotope_packet_count, self.packet_energy
                    )
                    * energy_factors
                )

                # rest frame gamma-ray energy and frequency
                # this probably works fine without the loop
                initial_packet_energies_rf = np.zeros(isotope_packet_count)
                initial_nus_rf = np.zeros(isotope_packet_count)
                for i in range(isotope_packet_count):
                    doppler_factor = doppler_factor_3d(
                        initial_directions[i],
                        initial_locations[i],
                        initial_times[i],
                    )
                    initial_packet_energies_rf[i] = (
                        initial_packet_energies_cmf[i] / doppler_factor
                    )
                    initial_nus_rf[i] = initial_nus_cmf[i] / doppler_factor

        locations = np.stack()

        return GXPacketCollection(
            locations,
            directions,
            packet_energies_rf,
            packet_energies_cmf,
            nus_rf,
            nus_cmf,
            shells,
            times,
        )

    @njit(**njit_dict_no_parallel)
    def calculate_positron_fraction(
        self, positron_energy, isotope_energy, isotope_intensity
    ):
        """Calculate the fraction of energy that an isotope
        releases as positron kinetic energy

        Parameters
        ----------
        positron_energy : float
            Average kinetic energy of positrons from decay
        isotope_energy : numpy array
            Photon energies released by the isotope
        isotope_intensity : numpy array
            Intensity of photon energy release

        Returns
        -------
        float
            Fraction of energy released as positron kinetic energy
        """
        return positron_energy / np.sum(isotope_energy * isotope_intensity)
