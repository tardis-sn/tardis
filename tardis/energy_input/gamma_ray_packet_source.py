import numpy as np
import pandas as pd

from tardis.energy_input.energy_source import (
    positronium_continuum,
)
from tardis.energy_input.GXPacket import (
    GXPacketCollection,
)
from tardis.energy_input.samplers import sample_energy
from tardis.energy_input.util import (
    H_CGS_KEV,
    doppler_factor_3d,
    get_index,
    get_random_unit_vector,
)
from tardis.transport.montecarlo.packet_source import BasePacketSource


class GammaRayPacketSource(BasePacketSource):
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
        self.average_power_per_mass = average_power_per_mass
        self.energy_plot_positron_rows = np.empty(0)
        super().__init__(**kwargs)

    def create_packet_mus(self, no_of_packets, *args, **kwargs):
        return super().create_packet_mus(no_of_packets, *args, **kwargs)

    def create_packet_radii(self, sampled_packets_df):
        """Initialize the random radii of packets in a shell

        Parameters
        ----------
        packet_count : int
            Number of packets in the shell
        sampled_packets_df : pd.DataFrame
            Dataframe where each row is a packet

        Returns
        -------
        array
            Array of length packet_count of random locations in the shell
        """
        z = np.random.random(len(sampled_packets_df))
        initial_radii = (
            z * sampled_packets_df["inner_velocity"] ** 3.0
            + (1.0 - z) * sampled_packets_df["outer_velocity"] ** 3.0
        ) ** (1.0 / 3.0)

        return initial_radii

    def create_packet_nus(
        self,
        no_of_packets,
        packets,
        positronium_fraction,
        positronium_energy,
        positronium_intensity,
    ):
        """Create an array of packet frequency-energies (i.e. E = h * nu)

        Parameters
        ----------
        no_of_packets : int
            Number of packets to produce frequency-energies for
        packets : pd.DataFrame
            DataFrame of packets
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
        """
        energy_array = np.zeros(no_of_packets)
        zs = np.random.random(no_of_packets)
        for i in range(no_of_packets):
            # positron
            if packets.iloc[i]["decay_type"] == "bp":
                # positronium formation 75% of the time if fraction is 1
                if zs[i] < positronium_fraction and np.random.random() < 0.75:
                    energy_array[i] = sample_energy(
                        positronium_energy, positronium_intensity
                    )
                else:
                    energy_array[i] = 511
            else:
                energy_array[i] = packets.iloc[i]["radiation_energy_kev"]

        return energy_array

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
        directions = np.zeros((3, no_of_packets))
        for i in range(no_of_packets):
            directions[:, i] = get_random_unit_vector()

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

    def create_packet_times_uniform_energy(
        self, no_of_packets, isotopes, decay_time
    ):
        """Samples the decay time from the mean lifetime of the isotopes

        Parameters
        ----------
        no_of_packets : int
            Number of packets
        isotopes : pd.Series
            Series of packet parent isotopes
        decay_time : array
            Series of packet decay time index

        Returns
        -------
        array
            Array of decay times
        """
        decay_times = np.zeros(len(no_of_packets))
        for i, isotope in enumerate(isotopes.to_numpy()):
            decay_time_min = self.times[decay_time[i]]
            if decay_time_min == self.times[-1]:
                decay_time_max = self.effective_times[-1]
            else:
                decay_time_max = self.times[decay_time[i] + 1]
            # rejection sampling
            while (decay_times[i] <= decay_time_min) or (
                decay_times[i] >= decay_time_max
            ):
                decay_times[i] = -self.taus[isotope] * np.log(
                    np.random.random()
                )
        return decay_times

    def create_packets(
        self, decays_per_isotope, number_of_packets, *args, **kwargs
    ):
        """Initialize a collection of GXPacket objects for the simulation
        to operate on.

        Parameters
        ----------
        decays_per_isotope : array int64
            Probability of decays per simulation shell per isotope per time step
        number_of_packets : int
            Number of packets to create

        Returns
        -------
        GXPacketCollection
        """
        # initialize arrays for most packet properties
        locations = np.zeros((3, number_of_packets))
        directions = np.zeros((3, number_of_packets))
        packet_energies_rf = np.zeros(number_of_packets)
        packet_energies_cmf = np.zeros(number_of_packets)
        nus_rf = np.zeros(number_of_packets)
        nus_cmf = np.zeros(number_of_packets)
        times = np.zeros(number_of_packets)
        # set packets to IN_PROCESS status
        statuses = np.ones(number_of_packets, dtype=np.int64) * 3

        self.energy_plot_positron_rows = np.zeros((number_of_packets, 4))

        # compute positronium continuum
        positronium_energy, positronium_intensity = positronium_continuum()

        # sample packets from dataframe, returning a dataframe where each row is
        # a sampled packet
        sampled_packets_df = decays_per_isotope.sample(
            n=number_of_packets,
            weights="decay_energy_erg",
            replace=True,
            random_state=np.random.RandomState(self.base_seed),
        )
        # get unique isotopes that have produced packets
        isotopes = pd.unique(sampled_packets_df.index.get_level_values(2))

        # compute the positron fraction for unique isotopes
        isotope_positron_fraction = self.calculate_positron_fraction(isotopes)

        # get the packet shell index
        shells = sampled_packets_df.index.get_level_values(1)

        # get the inner and outer velocity boundaries for each packet to compute
        # the initial radii
        sampled_packets_df["inner_velocity"] = self.inner_velocities[shells]
        sampled_packets_df["outer_velocity"] = self.outer_velocities[shells]

        # sample radii at time = 0
        initial_radii = self.create_packet_radii(sampled_packets_df)

        # get the time step index of the packets
        initial_time_indexes = sampled_packets_df.index.get_level_values(0)

        # get the time of the middle of the step for each packet
        packet_effective_times = np.array(
            [self.effective_times[i] for i in initial_time_indexes]
        )

        # packet decay time
        times = self.create_packet_times_uniform_energy(
            number_of_packets,
            sampled_packets_df.index.get_level_values(2),
            packet_effective_times,
        )

        # scale radius by packet decay time. This could be replaced with
        # Geometry object calculations. Note that this also adds a random
        # unit vector multiplication for 3D. May not be needed.
        locations = (
            initial_radii
            * packet_effective_times
            * self.create_packet_directions(number_of_packets)
        )

        # sample directions (valid at all times), non-relativistic
        directions = self.create_packet_directions(number_of_packets)

        # the individual gamma-ray energy that makes up a packet
        # co-moving frame, including positronium formation
        nu_energies_cmf = self.create_packet_nus(
            number_of_packets,
            sampled_packets_df,
            self.positronium_fraction,
            positronium_energy,
            positronium_intensity,
        )

        # equivalent frequencies
        nus_cmf = nu_energies_cmf / H_CGS_KEV

        # per packet co-moving frame total energy
        packet_energies_cmf = self.create_packet_energies(
            number_of_packets, self.packet_energy
        )

        # rest frame gamma-ray energy and frequency
        # this probably works fine without the loop
        # non-relativistic
        packet_energies_rf = np.zeros(number_of_packets)
        nus_rf = np.zeros(number_of_packets)
        for i in range(number_of_packets):
            doppler_factor = doppler_factor_3d(
                directions[:, i],
                locations[:, i],
                times[i],
            )
            packet_energies_rf[i] = packet_energies_cmf[i] / doppler_factor
            nus_rf[i] = nus_cmf[i] / doppler_factor

            # deposit positron energy in both output arrays
            # this is an average across all packets that are created
            # it could be changed to be only for packets that are from positrons
            self.energy_plot_positron_rows[i] = np.array(
                [
                    i,
                    isotope_positron_fraction[sampled_packets_df["isotopes"][i]]
                    * packet_energies_cmf[i],
                    # this needs to be sqrt(sum of squares) to get radius
                    np.linalg.norm(locations[i]),
                    times[i],
                ]
            )

            # this is an average across all packets that are created
            # it could be changed to be only for packets that are from positrons
            self.energy_df_rows[shells[i], times[i]] += (
                isotope_positron_fraction[sampled_packets_df["isotopes"][i]]
                * packet_energies_cmf[i]
            )

        return GXPacketCollection(
            locations,
            directions,
            packet_energies_rf,
            packet_energies_cmf,
            nus_rf,
            nus_cmf,
            statuses,
            shells,
            times,
        )

    def calculate_positron_fraction(self, isotopes):
        """Calculate the fraction of energy that an isotope
        releases as positron kinetic energy

        Parameters
        ----------
        isotopes : array
            Array of isotope names as strings

        Returns
        -------
        dict
            Fraction of energy released as positron kinetic energy per isotope
        """
        positron_fraction = {}

        for isotope in isotopes:
            isotope_energy = self.gamma_ray_lines[isotope][0, :]
            isotope_intensity = self.gamma_ray_lines[isotope][1, :]
            positron_fraction[isotope] = self.average_positron_energies[
                isotope
            ] / np.sum(isotope_energy * isotope_intensity)
        return positron_fraction
