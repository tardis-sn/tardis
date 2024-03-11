import numpy as np
from numba import njit

from tardis.energy_input.energy_source import (
    positronium_continuum,
)
from tardis.energy_input.GXPacket import (
    GXPacketCollection,
    initialize_packet_properties,
)
from tardis.energy_input.samplers import initial_packet_radius
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
        for k, shell in enumerate(decays_per_shell):
            initial_radii = initial_packet_radius(
                shell, self.inner_velocities[k], self.outer_velocities[k]
            )

            isotope_packet_count_df = decays_per_isotope.iloc[k]

            i = 0
            for (
                isotope_name,
                isotope_packet_count,
            ) in isotope_packet_count_df.items():
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

                for c in range(isotope_packet_count):
                    packet, decay_time_index = initialize_packet_properties(
                        isotope_energy,
                        isotope_intensity,
                        positronium_energy,
                        positronium_intensity,
                        self.positronium_fraction,
                        self.packet_energy,
                        k,
                        tau_start,
                        tau_end,
                        initial_radii[i],
                        self.times,
                        self.effective_times,
                        self.inventories[k],
                        self.average_power_per_mass,
                    )

                    self.energy_df_rows[k, decay_time_index] += (
                        isotope_positron_fraction * self.packet_energy * 1000
                    )

                    energy_plot_df_rows[packet_index] = np.array(
                        [
                            i,
                            packet.energy_rf,
                            packet.get_location_r(),
                            packet.time_current,
                            int(packet.status),
                            0,
                            0,
                            0,
                        ]
                    )

                    energy_plot_positron_rows[packet_index] = [
                        packet_index,
                        isotope_positron_fraction * self.packet_energy * 1000,
                        # * inv_volume_time[packet.shell, decay_time_index],
                        packet.get_location_r(),
                        packet.time_current,
                    ]

                    i += 1
                    packet_index += 1

        return (
            GXPacketCollection(),
            self.energy_df_rows,
            energy_plot_df_rows,
            energy_plot_positron_rows,
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
