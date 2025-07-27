import logging

import numpy as np
import pandas as pd

from tardis.energy_input.samplers import (
    PositroniumSampler,
)
from tardis.energy_input.transport.GXPacket import (
    GXPacketCollection,
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
        cumulative_decays_df: pd.DataFrame,
        isotope_decay_df: pd.DataFrame,
        positronium_fraction: float,
        inner_velocities: np.ndarray,
        outer_velocities: np.ndarray,
        times: np.ndarray,
        effective_times: np.ndarray,
        **kwargs,
    ) -> None:
        """
        Initialize gamma ray packet source.

        Initializes a gamma ray packet source for creating gamma ray packets
        from radioactive decay data, including support for positronium formation.

        Parameters
        ----------
        cumulative_decays_df : pd.DataFrame
            DataFrame containing cumulative decay data with columns including
            radiation type, decay energies, and multi-level indices for
            isotope, shell_number, and time_index.
        isotope_decay_df : pd.DataFrame
            DataFrame containing isotope decay data with decay constants
            and other isotope-specific parameters.
        positronium_fraction : float
            Fraction of positrons that form positronium (0.0 to 1.0).
            Used for modeling three-photon decay vs two-photon annihilation.
        inner_velocities : np.ndarray
            Array of inner shell velocities [cm/s] for each spatial shell.
        outer_velocities : np.ndarray
            Array of outer shell velocities [cm/s] for each spatial shell.
        times : np.ndarray
            Array of time steps [s] used in the simulation.
        effective_times : np.ndarray
            Array of effective time steps [s] accounting for simulation specifics.
        **kwargs
            Additional keyword arguments passed to the parent BasePacketSource class.

        Notes
        -----
        This packet source generates gamma ray packets from radioactive decay
        events, with proper handling of:
        - Spatial distribution within shells
        - Time-dependent decay processes
        - Positronium formation and decay modes
        - Doppler effects from moving material
        """
        self.cumulative_decays_df = cumulative_decays_df
        self.isotope_decay_df = isotope_decay_df
        self.positronium_fraction = positronium_fraction
        self.inner_velocities = inner_velocities
        self.outer_velocities = outer_velocities
        self.times = times
        self.effective_times = effective_times
        super().__init__(**kwargs)

    def create_packet_mus(self, no_of_packets: int, *args, **kwargs):
        """
        Create packet directional cosines.

        Creates packet directional cosines by calling the parent class method.
        This method is inherited from BasePacketSource.

        Parameters
        ----------
        no_of_packets : int
            Number of packets for which to create directional cosines.
        *args
            Variable length argument list passed to parent method.
        **kwargs
            Arbitrary keyword arguments passed to parent method.

        Returns
        -------
        The return value from the parent class create_packet_mus method.
        """
        return super().create_packet_mus(no_of_packets, *args, **kwargs)

    def create_packet_velocities(self, sampled_packets_df: pd.DataFrame) -> np.ndarray:
        """
        Initialize random radial velocities for packets within shells.

        Generates random initial velocities for packets distributed within
        spherical shells using a uniform distribution in volume.

        Parameters
        ----------
        sampled_packets_df : pd.DataFrame
            DataFrame where each row represents a packet, containing
            'inner_velocity' and 'outer_velocity' columns for shell boundaries.

        Returns
        -------
        np.ndarray
            Array of initial velocities [cm/s] with length equal to the number
            of packets in sampled_packets_df.

        Notes
        -----
        Uses the cube root method to ensure uniform distribution in volume:
        r^3 = z * r_inner^3 + (1-z) * r_outer^3, where z is uniform random [0,1].
        """
        np.random.seed(self.base_seed + 2 if self.base_seed is not None else None)
        z = np.random.random(len(sampled_packets_df))
        initial_velocities = (
            z * sampled_packets_df["inner_velocity"] ** 3.0
            + (1.0 - z) * sampled_packets_df["outer_velocity"] ** 3.0
        ) ** (1.0 / 3.0)

        return initial_velocities

    def create_packet_nus(
        self,
        packets: pd.DataFrame,
        positronium_fraction: float,
        number_of_packets: int,
    ) -> np.ndarray:
        """
        Create packet frequency-energies accounting for positronium formation.

        Generates an array of packet frequency-energies (E = h * nu) considering
        positronium formation and its decay modes for positron annihilation lines.

        Parameters
        ----------
        packets : pd.DataFrame
            DataFrame containing packet information with 'radiation_energy_keV' column.
        positronium_fraction : float
            Fraction of positrons that form positronium (0.0 to 1.0).
            Default is 0.0 for no positronium formation.
        number_of_packets : int
            Number of packets to generate frequency-energies for.

        Returns
        -------
        np.ndarray
            Array of sampled frequency-energies [keV] with length number_of_packets.

        Notes
        -----
        For positron annihilation lines (511 keV), this method:
        - Determines if positronium forms based on positronium_fraction
        - For ortho-positronium: samples from 3-photon decay spectrum
        - For para-positronium: uses the 511 keV line energy
        - For direct annihilation: uses the original 511 keV energy

        The para/ortho ratio is set by PARA_TO_ORTHO_RATIO constant (0.25).
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

    def create_packet_directions(
        self, no_of_packets: int, seed: int | None
    ) -> np.ndarray:
        """
        Create random isotropic directions for packets.

        Generates an array of random unit vectors representing isotropic
        directions for gamma ray packets.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to generate directions for.
        seed : int or None
            Random seed for reproducible direction generation.
            If None, uses current random state.

        Returns
        -------
        np.ndarray
            Array of shape (3, no_of_packets) containing unit direction vectors.
            Each column represents a 3D unit vector [x, y, z].

        Notes
        -----
        Directions are sampled uniformly on the unit sphere to ensure
        isotropic distribution in 3D space.
        """
        directions = get_random_unit_vectors(no_of_packets, seed)

        return directions

    def create_packet_energies(self, no_of_packets: int, energy: float) -> np.ndarray:
        """
        Create uniform packet energies for gamma ray packets.

        Generates an array of identical packet energies for a specified
        number of packets.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to create energies for.
        energy : float
            Energy value [erg] to assign to each packet.

        Returns
        -------
        np.ndarray
            Array of packet energies [erg] with length no_of_packets,
            where each element equals the input energy value.

        Notes
        -----
        This method creates uniform energy packets, where each packet
        carries the same energy regardless of the specific gamma ray
        line that created it. The total energy is conserved through
        the packet weighting system.
        """
        return np.ones(no_of_packets) * energy

    def create_packet_times_uniform_time(
        self, no_of_packets: int, start: float, end: float
    ) -> np.ndarray:
        """
        Sample packet decay times uniformly within a time interval.

        Generates decay times uniformly distributed between start and end times.
        This approach requires non-uniform packet energies to maintain
        energy conservation.

        Parameters
        ----------
        no_of_packets : int
            Number of packets to generate decay times for.
        start : float
            Start time [s] of the sampling interval.
        end : float
            End time [s] of the sampling interval.

        Returns
        -------
        np.ndarray
            Array of decay times [s] with length no_of_packets,
            uniformly distributed between start and end.

        Notes
        -----
        This method samples decay times uniformly in time, which means
        the packet energies must be weighted according to the decay
        rate at each time to properly represent the physical decay process.
        """
        z = np.random.random(no_of_packets)
        decay_times = z * start + (1 - z) * end
        return decay_times

    def create_packet_times_uniform_energy(
        self, no_of_packets: np.ndarray, isotopes: pd.Series, decay_time: np.ndarray
    ) -> np.ndarray:
        """
        Sample decay times from isotope mean lifetimes using rejection sampling.

        Generates decay times by sampling from exponential distributions
        based on isotope mean lifetimes, constrained to specific time intervals.

        Parameters
        ----------
        no_of_packets : np.ndarray
            Array indices for the packets (used for iteration).
        isotopes : pd.Series
            Series containing parent isotope names for each packet.
        decay_time : np.ndarray
            Array of time step indices indicating the time interval
            for each packet's decay.

        Returns
        -------
        np.ndarray
            Array of decay times [s] sampled from exponential distributions
            constrained to the appropriate time intervals.

        Notes
        -----
        This method uses rejection sampling to ensure decay times fall
        within the correct time bins. For each packet:
        1. Determines the time interval [t_min, t_max] from decay_time index
        2. Samples from exponential distribution: t = -tau * ln(random)
        3. Rejects and resamples if t is outside the interval

        Requires self.taus attribute containing isotope mean lifetimes.
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
        """
        Initialize a collection of gamma ray packets for simulation.

        Creates a collection of gamma ray packets from radioactive decay data,
        including proper spatial distribution, directional sampling, energy
        assignment, and Doppler corrections.

        Parameters
        ----------
        cumulative_decays_df : pd.DataFrame
            DataFrame containing cumulative decay data with columns including
            'radiation', 'decay_energy_erg', and multi-level index with
            'isotope', 'shell_number', and 'time_index'.
        number_of_packets : int
            Total number of gamma ray packets to create for the simulation.
        legacy_energy_per_packet : float, optional
            Legacy energy per packet [erg] for backwards compatibility.
            If None, energy per packet is calculated from total gamma ray
            energy divided by number of packets. Default is None.

        Returns
        -------
        GXPacketCollection
            Collection of gamma ray packets with initialized properties:
            - locations: 3D positions in the simulation domain
            - directions: isotropic unit direction vectors
            - energies: rest frame and comoving frame energies
            - frequencies: rest frame and comoving frame frequencies
            - metadata: shell numbers, decay times, source isotopes

        Notes
        -----
        The packet creation process includes:

        1. **Energy calculation**: Total gamma ray energy is divided equally
           among packets (uniform energy approach)
        2. **Spatial sampling**: Packets are distributed within shells based
           on decay energy weighting
        3. **Temporal placement**: Packets are positioned at decay times with
           appropriate radial expansion
        4. **Spectral sampling**: Frequencies include positronium formation
           effects for 511 keV annihilation lines
        5. **Doppler corrections**: Applied for relativistic motion between
           rest and comoving frames

        The method ensures energy conservation while providing proper
        statistical sampling of the decay process.
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


def legacy_calculate_positron_fraction(
    isotope_decay_df: pd.DataFrame, isotopes: np.ndarray, number_of_packets: int
) -> np.ndarray:
    """
    Calculate positron kinetic energy fraction relative to gamma ray energy.

    Computes the fraction of energy released as positron kinetic energy
    compared to gamma ray energy for each isotope associated with packets.

    Parameters
    ----------
    isotope_decay_df : pd.DataFrame
        DataFrame containing isotope decay data with multi-level index
        including 'isotope' and 'shell_number', and columns including
        'radiation', 'energy_per_channel_keV'.
    isotopes : np.ndarray
        Array of isotope names as strings, where each isotope is
        associated with a packet.
    number_of_packets : int
        Total number of gamma ray packets in the simulation.

    Returns
    -------
    np.ndarray
        Array of positron energy fractions with length number_of_packets.
        Each element represents the ratio of positron kinetic energy to
        gamma ray energy for the corresponding packet's source isotope.

    Notes
    -----
    This function:

    1. Filters decay data for shell_number == 0 to avoid double counting
    2. Separates gamma ray ('g') and beta plus ('bp') radiation channels
    3. Sums energy per isotope for each radiation type
    4. Calculates fraction = E_positron / E_gamma for each isotope
    5. Maps isotope fractions to packet array

    Isotopes not present in the decay DataFrame receive a fraction of 0.0.
    This is used for legacy compatibility and may be deprecated in favor
    of more sophisticated positron energy modeling.
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
