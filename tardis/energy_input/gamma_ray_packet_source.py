import logging

import numpy as np
import pandas as pd

from tardis.energy_input.transport.GXPacket import (
    GXPacketCollection,
)
from tardis.energy_input.samplers import (
    PositroniumSampler,
)
from tardis.energy_input.util import (
    H_CGS_KEV,
    doppler_factor_3D_all_packets,
    get_random_unit_vectors,
)
from tardis.transport.montecarlo.packet_source import BasePacketSource

logger = logging.getLogger(__name__)

POSITRON_ANNIHILATION_LINE = 511.0
PARA_TO_ORTHO_RATIO = 0.25


class GammaRayPacketSource(BasePacketSource):

    def __init__(
        self,
        cumulative_decays_df,
        isotope_decay_df,
        positronium_fraction,
        inner_velocities,
        outer_velocities,
        times,
        effective_times,
        **kwargs,
    ):
        """
        New Gamma ray packet source class

        Parameters
        ----------
        cumulative_decays_df : pd.DataFrame
            DataFrame containing the cumulative decay data
        isotope_decay_df : pd.DataFrame
            DataFrame of isotope decay data
        positronium_fraction : float
            Fraction of positrons that form positronium
        inner_velocities : array
            Array of inner shell velocities
        outer_velocities : array
            Array of outer shell velocities
        times : array
            Array of time steps
        effective_times : array
            Array of effective time steps
        """
        self.cumulative_decays_df = cumulative_decays_df
        self.isotope_decay_df = isotope_decay_df
        self.positronium_fraction = positronium_fraction
        self.inner_velocities = inner_velocities
        self.outer_velocities = outer_velocities
        self.times = times
        self.effective_times = effective_times
        super().__init__(**kwargs)

    def create_packet_mus(self, no_of_packets, *args, **kwargs):
        return super().create_packet_mus(no_of_packets, *args, **kwargs)

    def create_packet_velocities(self, sampled_packets_df):
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
        np.random.seed(self.base_seed + 2 if self.base_seed is not None else None)
        z = np.random.random(len(sampled_packets_df))
        initial_radii = (
            z * sampled_packets_df["inner_velocity"] ** 3.0
            + (1.0 - z) * sampled_packets_df["outer_velocity"] ** 3.0
        ) ** (1.0 / 3.0)

        return initial_radii

    def create_packet_nus(
        self,
        packets,
        positronium_fraction,
        number_of_packets,
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
            default is 0.0

        Returns
        -------
        array
            Array of sampled frequency-energies
        """
        energy_array = np.zeros(number_of_packets)

        all_packets = np.array([True] * number_of_packets)

        # positronium formation if fraction is greater than zero
        positronium_formation = (
            np.random.uniform(0, 1, number_of_packets) < positronium_fraction
        )
        # annihilation line of positrons
        annihilation_line = packets["radiation_energy_keV"] == POSITRON_ANNIHILATION_LINE
        # three photon decay of positronium
        three_photon_decay = np.random.random(number_of_packets) > PARA_TO_ORTHO_RATIO

        energy_array[all_packets] = packets.loc[
            all_packets, "radiation_energy_keV"
        ]

        energy_array[
            positronium_formation & annihilation_line & three_photon_decay
        ] = PositroniumSampler().sample_energy(
            samples=np.sum(
                positronium_formation & annihilation_line & three_photon_decay
            )
        )
        energy_array[
            positronium_formation & annihilation_line & ~three_photon_decay
        ] = POSITRON_ANNIHILATION_LINE

        return energy_array

    def create_packet_directions(self, no_of_packets, seed):
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
        directions = get_random_unit_vectors(no_of_packets, seed)

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
        self,
        cumulative_decays_df: pd.DataFrame,
        number_of_packets: int,
        legacy_energy_per_packet: float | None = None,
    ) -> GXPacketCollection:
        """Initialize a collection of GXPacket objects for the simulation.

        Parameters
        ----------
        cumulative_decays_df : pd.DataFrame
            DataFrame containing the cumulative decay data with columns including
            'radiation', 'decay_energy_erg', and multi-level index with 'isotope',
            'shell_number', and 'time_index'
        number_of_packets : int
            Number of gamma-ray packets to create
        legacy_energy_per_packet : float, optional
            Legacy energy per packet for backwards compatibility. If None,
            energy per packet is calculated from total gamma-ray energy
            divided by number of packets, by default None

        Returns
        -------
        GXPacketCollection
            Collection of gamma-ray packets with initialized properties including
            locations, directions, energies, frequencies, and metadata
        """
        # Calculate energy per packet
        if legacy_energy_per_packet is None:
            gamma_df = self.cumulative_decays_df[
                self.cumulative_decays_df["radiation"] == "g"
            ]
            total_energy_gamma = gamma_df["decay_energy_erg"].sum()
            energy_per_packet = total_energy_gamma / number_of_packets
        else:
            energy_per_packet = legacy_energy_per_packet
            total_energy_gamma = number_of_packets * legacy_energy_per_packet

        logger.info("Total energy in gamma-rays is %s", total_energy_gamma)
        logger.info("Energy per packet is %s", energy_per_packet)

        # initialize arrays for most packet properties
        locations = np.zeros((3, number_of_packets))
        directions = np.zeros((3, number_of_packets))
        packet_energies_rf = np.zeros(number_of_packets)
        packet_energies_cmf = np.zeros(number_of_packets)
        nus_rf = np.zeros(number_of_packets)
        nus_cmf = np.zeros(number_of_packets)
        statuses = np.ones(number_of_packets, dtype=np.int64) * 3

        # sample packets from the gamma-ray lines only (include X-rays!)
        sampled_packets_df_gamma = cumulative_decays_df[
            cumulative_decays_df["radiation"] == "g"
        ]

        # sample packets from the time evolving dataframe
        sampled_packets_df = sampled_packets_df_gamma.sample(
            n=number_of_packets,
            weights="decay_energy_erg",
            replace=True,
            random_state=np.random.RandomState(self.base_seed),
        )

        # get the isotopes and shells of the sampled packets
        source_isotopes = sampled_packets_df.index.get_level_values("isotope")
        shells = sampled_packets_df.index.get_level_values("shell_number")

        # get the inner and outer velocity boundaries for each packet to compute
        sampled_packets_df["inner_velocity"] = self.inner_velocities[shells]
        sampled_packets_df["outer_velocity"] = self.outer_velocities[shells]

        # The radii of the packets at what ever time they are emitted
        initial_velocities = self.create_packet_velocities(sampled_packets_df)

        # get the time step index of the packets
        decay_time_indices = sampled_packets_df.index.get_level_values("time_index")

        effective_decay_times = self.times[decay_time_indices]

        # scale radius by packet decay time. This could be replaced with
        # Geometry object calculations. Note that this also adds a random
        # unit vector multiplication for 3D. May not be needed.
        locations = (
            initial_velocities.values
            * effective_decay_times
            * self.create_packet_directions(number_of_packets, seed=self.base_seed)
        )

        # sample directions (valid at all times), non-relativistic
        # the seed is changed to not have packets that are all going outwards as the
        # create_packet_directions method is also used for the location sampling
        directions_seed = self.base_seed + 1 if self.base_seed is not None else None
        directions = self.create_packet_directions(
            number_of_packets, seed=directions_seed
        )

        # the individual gamma-ray energy that makes up a packet
        # co-moving frame, including positronium formation
        nu_energies_cmf = self.create_packet_nus(
            sampled_packets_df,
            self.positronium_fraction,
            number_of_packets,
        )

        nus_cmf = nu_energies_cmf / H_CGS_KEV

        packet_energies_cmf = self.create_packet_energies(
            number_of_packets, energy_per_packet
        )
        packet_energies_rf = np.zeros(number_of_packets)
        nus_rf = np.zeros(number_of_packets)

        doppler_factors = doppler_factor_3D_all_packets(
            directions, locations, effective_decay_times
        )

        packet_energies_rf = packet_energies_cmf / doppler_factors
        nus_rf = nus_cmf / doppler_factors

        return GXPacketCollection(
            locations,
            directions,
            packet_energies_rf,
            packet_energies_cmf,
            nus_rf,
            nus_cmf,
            statuses,
            shells,
            effective_decay_times,
            decay_time_indices,
            source_isotopes=source_isotopes,
        )


def legacy_calculate_positron_fraction(isotope_decay_df, isotopes, number_of_packets):
    """Calculate the fraction of energy that an isotope
    releases as positron kinetic energy compared to gamma-ray energy

    Parameters
    ----------
    isotope_decay_df : pd.DataFrame
        DataFrame of isotope decay data
    isotopes : array
        Array of isotope names as strings. Here each isotope is associated with a packet.
    number_of_packets : int
        Number of gamma-ray packets

    Returns
    -------
    dict
        Fraction of energy released as positron kinetic energy per isotope
    """
    isotope_positron_fraction = np.zeros(number_of_packets)

    # Find the positron fraction from the zeroth shell of the dataframe
    # this is because the total positron kinetic energy is the same for all shells
    shell_number_0 = isotope_decay_df[
        isotope_decay_df.index.get_level_values("shell_number") == 0
    ]

    gamma_decay_df = shell_number_0[shell_number_0["radiation"] == "g"]

    positrons_decay_df = shell_number_0[shell_number_0["radiation"] == "bp"]
    # Find the total energy released from positrons per isotope from the dataframe
    positron_energy_per_isotope = positrons_decay_df.groupby("isotope")[
        "energy_per_channel_keV"
    ].sum()
    # Find the total energy released from gamma-ray per isotope from the dataframe
    # TODO: Can be tested with total energy released from all radiation types
    gamma_energy_per_isotope = gamma_decay_df.groupby("isotope")[
        "energy_per_channel_keV"
    ].sum()
    # TODO: Possibly move this for loop
    for i, isotope in enumerate(isotopes):
        if (
            isotope in positron_energy_per_isotope
        ):  # check if isotope is in the dataframe
            isotope_positron_fraction[i] = (
                positron_energy_per_isotope[isotope] / gamma_energy_per_isotope[isotope]
            )
    return isotope_positron_fraction
