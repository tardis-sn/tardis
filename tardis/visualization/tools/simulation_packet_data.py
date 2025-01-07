import astropy.units as u
import pandas as pd


class SimulationPacketData:
    """
    The data structure representing simulation packet properties.

    This class preprocesses simulation packet data, encapsulating details about packet
    interactions, emission, and absorption for analysis and visualization.
    """

    def __init__(
        self,
        last_interaction_type,
        last_line_interaction_in_id,
        last_line_interaction_out_id,
        last_line_interaction_in_nu,
        last_interaction_in_r,
        lines_df,
        packet_nus,
        packet_energies,
        r_inner,
        spectrum_delta_frequency,
        spectrum_frequency_bins,  # stores _frequency (bin edges) not frequency
        spectrum_luminosity_density_lambda,
        spectrum_wavelength,
        t_inner,
        time_of_simulation,
    ):
        """
        Initialize the SimulationPacketData with required properties of simulation model.

        Parameters
        ----------
        last_interaction_type : np.array
            Interaction type (no-interaction: -1, e-scattering: 1 and
            line interaction: 2) values of emitted packets
        last_line_interaction_in_id : np.array
            IDs of atomic lines with which emitted packet had their last
            absorption (interaction in)
        last_line_interaction_out_id : np.array
            IDs of atomic lines with which emitted packet had their last
            emission (interaction out)
        last_line_interaction_in_nu : np.array
            Frequency values of the last absorption of emitted packets
        last_line_interaction_in_r : np.array
            Radius of the last interaction experienced by emitted packets
        lines_df : pd.DataFrame
            Data about the atomic lines present in simulation model's plasma
        packet_nus : astropy.Quantity
            Frequency values of the last emission of emitted packets, having
            unit of Hz
        packet_energies : astropy.Quantity
            Energy values of emitted packets, having unit of erg
        r_inner : astropy.Quantity
            Radius of innermost shell, having unit of cm
        spectrum_delta_frequency : astropy.Quantity
            Frequency bin width of spectrum, having unit of Hz
        spectrum_frequency_bins : astropy.Quantity
            Frequency bin edges of spectrum, having unit of Hz
        spectrum_wavelength : astropy.Quantity
            Wavelength values of spectrum, having unit of Angstrom
        t_inner : astropy.Quantity
            Temperature of innermost shell, having unit of K
        time_of_simulation : astropy.Quantity
            Time of simulation, having unit of s (second)
        """
        # Save packets properties in a dataframe for easier data manipulation
        packet_nus = u.Quantity(packet_nus, u.Hz)
        self.packets_df = pd.DataFrame(
            {
                "nus": packet_nus,
                "lambdas": packet_nus.to("angstrom", u.spectral()),
                "energies": packet_energies,
                "last_interaction_type": last_interaction_type,
                "last_line_interaction_out_id": last_line_interaction_out_id,
                "last_line_interaction_in_id": last_line_interaction_in_id,
                "last_line_interaction_in_nu": last_line_interaction_in_nu,
                "last_interaction_in_r": last_interaction_in_r,
            }
        )

        # Save other properties
        self.lines_df = lines_df
        self.r_inner = r_inner
        self.spectrum_delta_frequency = spectrum_delta_frequency
        self.spectrum_frequency_bins = spectrum_frequency_bins
        self.spectrum_frequency = spectrum_frequency_bins[:-1]
        self.spectrum_luminosity_density_lambda = (
            spectrum_luminosity_density_lambda
        )
        self.spectrum_wavelength = spectrum_wavelength
        self.t_inner = t_inner
        self.time_of_simulation = time_of_simulation

        # Create dataframe of packets that experience line interaction
        line_mask = (self.packets_df["last_interaction_type"] > -1) & (
            self.packets_df["last_line_interaction_in_id"] > -1
        )  # & operator is quite faster than np.logical_and on pd.Series
        self.packets_df_line_interaction = self.packets_df.loc[line_mask].copy()

        # Add columns for atomic number of last interaction out
        self.packets_df_line_interaction["last_line_interaction_atom"] = (
            self.lines_df["atomic_number"]
            .iloc[
                self.packets_df_line_interaction["last_line_interaction_out_id"]
            ]
            .to_numpy()
        )
        # Add columns for the species id of last interaction
        # Species id is given by 100 * Z + X, where Z is atomic number and X is ion number
        self.packets_df_line_interaction["last_line_interaction_species"] = (
            self.lines_df["atomic_number"]
            .iloc[
                self.packets_df_line_interaction["last_line_interaction_out_id"]
            ]
            .to_numpy()
            * 100
            + self.lines_df["ion_number"]
            .iloc[
                self.packets_df_line_interaction["last_line_interaction_out_id"]
            ]
            .to_numpy()
        )

    @classmethod
    def from_simulation(cls, sim, packets_mode):
        """
        Create an instance of SimulationPacketData from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation
        packets_mode : {'virtual', 'real'}
            Mode of packets to be considered, either real or virtual

        Returns
        -------
        SimulationPacketData
        """
        if packets_mode not in ["virtual", "real"]:
            raise ValueError(
                "Invalid value passed to packets_mode. Only "
                "allowed values are 'virtual' or 'real'"
            )
        # Properties common among both packet modes
        lines_df = sim.plasma.atomic_data.lines.reset_index().set_index(
            "line_id"
        )
        transport_state = sim.transport.transport_state
        r_inner = sim.simulation_state.geometry.r_inner_active
        t_inner = sim.simulation_state.packet_source.temperature
        time_of_simulation = (
            transport_state.packet_collection.time_of_simulation * u.s
        )
        spectrum = (
            sim.spectrum_solver.spectrum_virtual_packets
            if packets_mode == "virtual"
            else sim.spectrum_solver.spectrum_real_packets
        )

        if packets_mode == "virtual":
            vpacket_tracker = transport_state.vpacket_tracker
            return cls(
                last_interaction_type=vpacket_tracker.last_interaction_type,
                last_line_interaction_in_id=vpacket_tracker.last_interaction_in_id,
                last_line_interaction_out_id=vpacket_tracker.last_interaction_out_id,
                last_line_interaction_in_nu=vpacket_tracker.last_interaction_in_nu,
                last_interaction_in_r=vpacket_tracker.last_interaction_in_r,
                lines_df=lines_df,
                packet_nus=u.Quantity(vpacket_tracker.nus, "Hz"),
                packet_energies=u.Quantity(vpacket_tracker.energies, "erg"),
                r_inner=r_inner,
                spectrum_delta_frequency=spectrum.delta_frequency,
                spectrum_frequency_bins=spectrum._frequency,
                spectrum_luminosity_density_lambda=spectrum.luminosity_density_lambda,
                spectrum_wavelength=spectrum.wavelength,
                t_inner=t_inner,
                time_of_simulation=time_of_simulation,
            )
        else: # real packets
            # Packets-specific properties need to be only for those packets
            # which got emitted
            mask = transport_state.emitted_packet_mask
            return cls(
                last_interaction_type=transport_state.last_interaction_type[mask],
                last_line_interaction_in_id=transport_state.last_line_interaction_in_id[
                    mask
                ],
                last_line_interaction_out_id=transport_state.last_line_interaction_out_id[
                    mask
                ],
                last_line_interaction_in_nu=transport_state.last_interaction_in_nu[
                    mask
                ],
                last_interaction_in_r=transport_state.last_interaction_in_r[mask],
                lines_df=lines_df,
                packet_nus=transport_state.packet_collection.output_nus[mask],
                packet_energies=transport_state.packet_collection.output_energies[
                    mask
                ],
                r_inner=r_inner,
                spectrum_delta_frequency=spectrum.delta_frequency,
                spectrum_frequency_bins=spectrum._frequency,
                spectrum_luminosity_density_lambda=spectrum.luminosity_density_lambda,
                spectrum_wavelength=spectrum.wavelength,
                t_inner=t_inner,
                time_of_simulation=time_of_simulation,
            )

    @classmethod
    def from_hdf(cls, hdf_fpath, packets_mode):
        """
        Create an instance of SimulationPacketData from a simulation HDF file.

        Parameters
        ----------
        hdf_fpath : str
            Valid path to the HDF file where simulation is saved
        packets_mode : {'virtual', 'real'}
            Mode of packets to be considered, either real or virtual

        Returns
        -------
        SimulationPacketData
        """
        if packets_mode not in ["virtual", "real"]:
            raise ValueError(
                "Invalid value passed to packets_mode. Only "
                "allowed values are 'virtual' or 'real'"
            )
        with pd.HDFStore(hdf_fpath, "r") as hdf:
            lines_df = (
                hdf["/simulation/plasma/lines"]
                .reset_index()
                .set_index("line_id")
            )
            r_inner = u.Quantity(
                hdf["/simulation/simulation_state/r_inner"].to_numpy(), "cm"
            )  # Convert pd.Series to np.array to construct quantity from it
            t_inner = u.Quantity(
                hdf["/simulation/simulation_state/scalars"].t_inner, "K"
            )
            time_of_simulation = u.Quantity(
                hdf[
                    "/simulation/transport/transport_state/scalars"
                ].time_of_simulation,
                "s",
            )

        spectrum_prefix = (
            f"/simulation/spectrum_solver/spectrum_{packets_mode}_packets"
        )

        if packets_mode == "virtual":
            packet_prefix = "/simulation/transport/transport_state/virt_packet"
            return cls(
                last_interaction_type=hdf[
                    f"{packet_prefix}_last_interaction_type"
                ],
                last_line_interaction_in_id=hdf[
                    f"{packet_prefix}_last_line_interaction_in_id"
                ],
                last_line_interaction_out_id=hdf[
                    f"{packet_prefix}_last_line_interaction_out_id"
                ],
                last_line_interaction_in_nu=u.Quantity(
                    hdf[f"{packet_prefix}_last_interaction_in_nu"].to_numpy(),
                    "Hz",
                ),
                last_interaction_in_r=u.Quantity(
                    hdf[f"{packet_prefix}_last_interaction_in_r"].to_numpy(),
                    "cm",
                ),
                lines_df=lines_df,
                packet_nus=u.Quantity(
                    hdf[f"{packet_prefix}_nus"].to_numpy(), "Hz"
                ),
                packet_energies=u.Quantity(
                    hdf[f"{packet_prefix}_energies"].to_numpy(), "erg"
                ),
                r_inner=r_inner,
                spectrum_delta_frequency=u.Quantity(
                    hdf[f"{spectrum_prefix}/scalars"].delta_frequency, "Hz"
                ),
                spectrum_frequency_bins=u.Quantity(
                    hdf[f"{spectrum_prefix}/_frequency"].to_numpy(), "Hz"
                ),
                spectrum_luminosity_density_lambda=u.Quantity(
                    hdf[
                        f"{spectrum_prefix}/luminosity_density_lambda"
                    ].to_numpy(),
                    "erg / s cm",
                ).to("erg / s AA"),
                spectrum_wavelength=u.Quantity(
                    hdf[f"{spectrum_prefix}/wavelength"].to_numpy(), "cm"
                ).to("AA"),
                t_inner=t_inner,
                time_of_simulation=time_of_simulation,
            )
        else: # real packets
            emitted_packet_mask = hdf[
                "/simulation/transport/transport_state/emitted_packet_mask"
            ].to_numpy()
            packet_prefix = "/simulation/transport/transport_state"
            return cls(
                last_interaction_type=hdf[
                    f"{packet_prefix}/last_interaction_type"
                ].to_numpy()[emitted_packet_mask],
                last_line_interaction_in_id=hdf[
                    f"{packet_prefix}/last_line_interaction_in_id"
                ].to_numpy()[emitted_packet_mask],
                last_line_interaction_out_id=hdf[
                    f"{packet_prefix}/last_line_interaction_out_id"
                ].to_numpy()[emitted_packet_mask],
                last_line_interaction_in_nu=u.Quantity(
                    hdf[f"{packet_prefix}/last_interaction_in_nu"].to_numpy()[
                        emitted_packet_mask
                    ],
                    "Hz",
                ),
                last_interaction_in_r=u.Quantity(
                    hdf[f"{packet_prefix}/last_interaction_in_r"].to_numpy()[
                        emitted_packet_mask
                    ],
                    "cm",
                ),
                lines_df=lines_df,
                packet_nus=u.Quantity(
                    hdf[f"{packet_prefix}/output_nu"].to_numpy()[
                        emitted_packet_mask
                    ],
                    "Hz",
                ),
                packet_energies=u.Quantity(
                    hdf[f"{packet_prefix}/output_energy"].to_numpy()[
                        emitted_packet_mask
                    ],
                    "erg",
                ),
                r_inner=r_inner,
                spectrum_delta_frequency=u.Quantity(
                    hdf[f"{spectrum_prefix}/scalars"].delta_frequency, "Hz"
                ),
                spectrum_frequency_bins=u.Quantity(
                    hdf[f"{spectrum_prefix}/_frequency"].to_numpy(), "Hz"
                ),
                spectrum_luminosity_density_lambda=u.Quantity(
                    hdf[f"{spectrum_prefix}/luminosity_density_lambda"].to_numpy(),
                    "erg / s cm",
                ).to("erg / s AA"),
                spectrum_wavelength=u.Quantity(
                    hdf[f"{spectrum_prefix}/wavelength"].to_numpy(), "cm"
                ).to("AA"),
                t_inner=t_inner,
                time_of_simulation=time_of_simulation,
            )
