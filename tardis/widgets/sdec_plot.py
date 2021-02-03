"""
Spectral element DEComposition (SDEC) Plot for TARDIS simulation models.

This plot is a spectral diagnostics plot similar to those originally
proposed by M. Kromer (see, for example, Kromer et al. 2013, figure 4).
"""
import numpy as np
import pandas as pd
import astropy.units as u
import astropy.modeling.blackbody as abb

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as clr
import plotly.graph_objects as go

from tardis.util.base import atomic_number2element_symbol
from tardis.widgets import plot_util as pu


class SDECData:
    """The data of simulation model used by Spectral element DEComposition (SDEC) Plot.

    This preprocesses the data required by SDECPlotter class for doing
    calculations and plotting.
    """

    def __init__(
        self,
        last_interaction_type,
        last_line_interaction_in_id,
        last_line_interaction_out_id,
        last_line_interaction_in_nu,
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
        Initialize the SDECData with required properties of simulation model.

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
        self.packets_df = pd.DataFrame(
            {
                "nus": packet_nus,
                "lambdas": packet_nus.to("angstrom", u.spectral()),
                "energies": packet_energies,
                "last_interaction_type": last_interaction_type,
                "last_line_interaction_out_id": last_line_interaction_out_id,
                "last_line_interaction_in_id": last_line_interaction_in_id,
                "last_line_interaction_in_nu": last_line_interaction_in_nu,
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

        # Add columns for atomic number of last interaction in/out
        self.packets_df_line_interaction["last_line_interaction_out_atom"] = (
            self.lines_df["atomic_number"]
            .iloc[
                self.packets_df_line_interaction["last_line_interaction_out_id"]
            ]
            .to_numpy()
        )
        self.packets_df_line_interaction["last_line_interaction_in_atom"] = (
            self.lines_df["atomic_number"]
            .iloc[
                self.packets_df_line_interaction["last_line_interaction_in_id"]
            ]
            .to_numpy()
        )

    @classmethod
    def from_simulation(cls, sim, packets_mode):
        """
        Create an instance of SDECData from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation
        packets_mode : {'virtual', 'real'}
            Mode of packets to be considered, either real or virtual

        Returns
        -------
        SDECData
        """
        # Properties common among both packet modes
        lines_df = sim.plasma.atomic_data.lines.reset_index().set_index(
            "line_id"
        )
        r_inner = sim.model.r_inner
        t_inner = sim.model.t_inner
        time_of_simulation = sim.runner.time_of_simulation

        if packets_mode == "virtual":
            return cls(
                last_interaction_type=sim.runner.virt_packet_last_interaction_type,
                last_line_interaction_in_id=sim.runner.virt_packet_last_line_interaction_in_id,
                last_line_interaction_out_id=sim.runner.virt_packet_last_line_interaction_out_id,
                last_line_interaction_in_nu=sim.runner.virt_packet_last_interaction_in_nu,
                lines_df=lines_df,
                packet_nus=u.Quantity(sim.runner.virt_packet_nus, "Hz"),
                packet_energies=u.Quantity(
                    sim.runner.virt_packet_energies, "erg"
                ),
                r_inner=r_inner,
                spectrum_delta_frequency=sim.runner.spectrum_virtual.delta_frequency,
                spectrum_frequency_bins=sim.runner.spectrum_virtual._frequency,
                spectrum_luminosity_density_lambda=sim.runner.spectrum_virtual.luminosity_density_lambda,
                spectrum_wavelength=sim.runner.spectrum_virtual.wavelength,
                t_inner=t_inner,
                time_of_simulation=time_of_simulation,
            )

        elif packets_mode == "real":
            # Packets-specific properties need to be only for those packets
            # which got emitted
            return cls(
                last_interaction_type=sim.runner.last_interaction_type[
                    sim.runner.emitted_packet_mask
                ],
                last_line_interaction_in_id=sim.runner.last_line_interaction_in_id[
                    sim.runner.emitted_packet_mask
                ],
                last_line_interaction_out_id=sim.runner.last_line_interaction_out_id[
                    sim.runner.emitted_packet_mask
                ],
                last_line_interaction_in_nu=sim.runner.last_interaction_in_nu[
                    sim.runner.emitted_packet_mask
                ],
                lines_df=lines_df,
                packet_nus=sim.runner.output_nu[sim.runner.emitted_packet_mask],
                packet_energies=sim.runner.output_energy[
                    sim.runner.emitted_packet_mask
                ],
                r_inner=r_inner,
                spectrum_delta_frequency=sim.runner.spectrum.delta_frequency,
                spectrum_frequency_bins=sim.runner.spectrum._frequency,
                spectrum_luminosity_density_lambda=sim.runner.spectrum.luminosity_density_lambda,
                spectrum_wavelength=sim.runner.spectrum.wavelength,
                t_inner=t_inner,
                time_of_simulation=time_of_simulation,
            )
        else:
            raise ValueError(
                "Invalid value passed to packets_mode. Only "
                "allowed values are 'virtual' or 'real'"
            )

    @classmethod
    def from_hdf(cls, hdf_fpath, packets_mode):
        """
        Create an instance of SDECData from a simulation HDF file.

        Parameters
        ----------
        hdf_fpath : str
            Valid path to the HDF file where simulation is saved
        packets_mode : {'virtual', 'real'}
            Mode of packets to be considered, either real or virtual

        Returns
        -------
        SDECData
        """
        with pd.HDFStore(hdf_fpath, "r") as hdf:
            lines_df = (
                hdf["/simulation/plasma/lines"]
                .reset_index()
                .set_index("line_id")
            )
            r_inner = u.Quantity(
                hdf["/simulation/model/r_inner"].to_numpy(), "cm"
            )  # Convert pd.Series to np.array to construct quantity from it
            t_inner = u.Quantity(hdf["/simulation/model/scalars"].t_inner, "K")
            time_of_simulation = u.Quantity(
                hdf["/simulation/runner/scalars"].time_of_simulation, "s"
            )

            if packets_mode == "virtual":
                return cls(
                    last_interaction_type=hdf[
                        "/simulation/runner/virt_packet_last_interaction_type"
                    ],
                    last_line_interaction_in_id=hdf[
                        "/simulation/runner/virt_packet_last_line_interaction_in_id"
                    ],
                    last_line_interaction_out_id=hdf[
                        "/simulation/runner/virt_packet_last_line_interaction_out_id"
                    ],
                    last_line_interaction_in_nu=u.Quantity(
                        hdf[
                            "/simulation/runner/virt_packet_last_interaction_in_nu"
                        ].to_numpy(),
                        "Hz",
                    ),
                    lines_df=lines_df,
                    packet_nus=u.Quantity(
                        hdf["/simulation/runner/virt_packet_nus"].to_numpy(),
                        "Hz",
                    ),
                    packet_energies=u.Quantity(
                        hdf[
                            "/simulation/runner/virt_packet_energies"
                        ].to_numpy(),
                        "erg",
                    ),
                    r_inner=r_inner,
                    spectrum_delta_frequency=u.Quantity(
                        hdf[
                            "/simulation/runner/spectrum_virtual/scalars"
                        ].delta_frequency,
                        "Hz",
                    ),
                    spectrum_frequency_bins=u.Quantity(
                        hdf[
                            "/simulation/runner/spectrum_virtual/_frequency"
                        ].to_numpy(),
                        "Hz",
                    ),
                    spectrum_luminosity_density_lambda=u.Quantity(
                        hdf[
                            "/simulation/runner/spectrum_virtual/luminosity_density_lambda"
                        ].to_numpy(),
                        "erg / s cm",  # luminosity_density_lambda is saved in hdf in CGS
                    ).to("erg / s AA"),
                    spectrum_wavelength=u.Quantity(
                        hdf[
                            "/simulation/runner/spectrum_virtual/wavelength"
                        ].to_numpy(),
                        "cm",  # wavelength is saved in hdf in CGS
                    ).to("AA"),
                    t_inner=t_inner,
                    time_of_simulation=time_of_simulation,
                )

            elif packets_mode == "real":
                emitted_packet_mask = hdf[
                    "/simulation/runner/emitted_packet_mask"
                ].to_numpy()
                return cls(
                    # First convert series read from hdf to array before masking
                    # to eliminate index info which creates problems otherwise
                    last_interaction_type=hdf[
                        "/simulation/runner/last_interaction_type"
                    ].to_numpy()[emitted_packet_mask],
                    last_line_interaction_in_id=hdf[
                        "/simulation/runner/last_line_interaction_in_id"
                    ].to_numpy()[emitted_packet_mask],
                    last_line_interaction_out_id=hdf[
                        "/simulation/runner/last_line_interaction_out_id"
                    ].to_numpy()[emitted_packet_mask],
                    last_line_interaction_in_nu=u.Quantity(
                        hdf[
                            "/simulation/runner/last_interaction_in_nu"
                        ].to_numpy()[emitted_packet_mask],
                        "Hz",
                    ),
                    lines_df=lines_df,
                    packet_nus=u.Quantity(
                        hdf["/simulation/runner/output_nu"].to_numpy()[
                            emitted_packet_mask
                        ],
                        "Hz",
                    ),
                    packet_energies=u.Quantity(
                        hdf["/simulation/runner/output_energy"].to_numpy()[
                            emitted_packet_mask
                        ],
                        "erg",
                    ),
                    r_inner=r_inner,
                    spectrum_delta_frequency=u.Quantity(
                        hdf[
                            "/simulation/runner/spectrum/scalars"
                        ].delta_frequency,
                        "Hz",
                    ),
                    spectrum_frequency_bins=u.Quantity(
                        hdf[
                            "/simulation/runner/spectrum/_frequency"
                        ].to_numpy(),
                        "Hz",
                    ),
                    spectrum_luminosity_density_lambda=u.Quantity(
                        hdf[
                            "/simulation/runner/spectrum/luminosity_density_lambda"
                        ].to_numpy(),
                        "erg / s cm",
                    ).to("erg / s AA"),
                    spectrum_wavelength=u.Quantity(
                        hdf[
                            "/simulation/runner/spectrum/wavelength"
                        ].to_numpy(),
                        "cm",
                    ).to("AA"),
                    t_inner=t_inner,
                    time_of_simulation=time_of_simulation,
                )
            else:
                raise ValueError(
                    "Invalid value passed to packets_mode. Only "
                    "allowed values are 'virtual' or 'real'"
                )


class SDECPlotter:
    """
    Plotting interface for Spectral element DEComposition (SDEC) Plot.

    It performs necessary calculations to generate SDEC Plot for a simulation
    model, and allows to plot it in matplotlib and plotly.
    """

    def __init__(self, data):
        """
        Initialize the SDECPlotter with required data of simulation model.

        Parameters
        ----------
        data : dict of SDECData
            Dictionary to store data required for SDEC plot, for both packet
            modes i.e. real and virtual
        """
        self.data = data

    @classmethod
    def from_simulation(cls, sim):
        """
        Create an instance of SDECPlotter from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation

        Returns
        -------
        SDECPlotter
        """
        return cls(
            dict(
                virtual=SDECData.from_simulation(sim, "virtual"),
                real=SDECData.from_simulation(sim, "real"),
            )
        )

    @classmethod
    def from_hdf(cls, hdf_fpath):
        """
        Create an instance of SDECPlotter from a simulation HDF file.

        Parameters
        ----------
        hdf_fpath : str
            Valid path to the HDF file where simulation is saved

        Returns
        -------
        SDECPlotter
        """
        return cls(
            dict(
                virtual=SDECData.from_hdf(hdf_fpath, "virtual"),
                real=SDECData.from_hdf(hdf_fpath, "real"),
            )
        )

    def _calculate_plotting_data(
        self, packets_mode, packet_wvl_range, distance
    ):
        """
        Calculate data to be used in plotting based on parameters passed.

        Parameters
        ----------
        packets_mode : {'virtual', 'real'}
            Mode of packets to be considered, either real or virtual
        packet_wvl_range : astropy.Quantity
            Wavelength range to restrict the analysis of escaped packets. It
            should be a quantity having units of Angstrom, containing two
            values - lower lambda and upper lambda i.e.
            [lower_lambda, upper_lambda] * u.AA
        distance : astropy.Quantity
            Distance used to calculate flux instead of luminosities in the plot.
            It should have a length unit like m, Mpc, etc.

        Notes
        -----
        It doesn't return the calculated properties but save them in instance
        itself. So it should be always called before starting plotting to
        update the plotting data based on parameters passed.
        """
        if packets_mode not in ["virtual", "real"]:
            raise ValueError(
                "Invalid value passed to packets_mode. Only "
                "allowed values are 'virtual' or 'real'"
            )

        # Store the plottable range of each spectrum property which is
        # same as entire range, initially
        self.plot_frequency_bins = self.data[
            packets_mode
        ].spectrum_frequency_bins
        self.plot_wavelength = self.data[packets_mode].spectrum_wavelength
        self.plot_frequency = self.data[packets_mode].spectrum_frequency

        # Filter their plottable range based on packet_wvl_range specified
        if packet_wvl_range:
            packet_nu_range = packet_wvl_range.to("Hz", u.spectral())

            # Index of value just before the 1st value that is > packet_nu_range[1]
            start_idx = (
                np.argmax(self.plot_frequency_bins > packet_nu_range[1]) - 1
            )
            # Index of value just after the last value that is < packet_nu_range[0]
            end_idx = np.argmin(self.plot_frequency_bins < packet_nu_range[0])
            self.plot_frequency_bins = self.plot_frequency_bins[
                start_idx : end_idx + 1
            ]

            # Since spectrum frequency (& hence wavelength) were created from
            # frequency_bins[:-1], so we exclude end_idx when creating the mask
            self.packet_wvl_range_mask = np.zeros(
                self.plot_wavelength.size, dtype=bool
            )
            self.packet_wvl_range_mask[start_idx:end_idx] = True

            self.plot_wavelength = self.plot_wavelength[
                self.packet_wvl_range_mask
            ]
            self.plot_frequency = self.plot_frequency[
                self.packet_wvl_range_mask
            ]
        else:
            self.packet_wvl_range_mask = np.ones(
                self.plot_wavelength.size, dtype=bool
            )

        # Make sure number of bin edges are always one more than wavelengths
        assert self.plot_frequency_bins.size == self.plot_wavelength.size + 1

        # Calculate the area term to convert luminosity to flux
        if distance is None:
            self.lum_to_flux = 1  # so that this term will have no effect
        else:
            if distance <= 0:
                raise ValueError(
                    "distance passed must be greater than 0. If you intended "
                    "to plot luminosities instead of flux, set distance=None "
                    "or don't specify distance parameter in the function call."
                )
            else:
                self.lum_to_flux = 4.0 * np.pi * (distance.to("cm")) ** 2

        # Calculate luminosities to be shown in plot
        (
            self.emission_luminosities_df,
            self.elements,
        ) = self._calculate_emission_luminosities(
            packets_mode=packets_mode, packet_wvl_range=packet_wvl_range
        )
        self.absorption_luminosities_df = (
            self._calculate_absorption_luminosities(
                packets_mode=packets_mode, packet_wvl_range=packet_wvl_range
            )
        )
        self.photosphere_luminosity = self._calculate_photosphere_luminosity(
            packets_mode=packets_mode
        )
        self.modeled_spectrum_luminosity = (
            self.data[packets_mode].spectrum_luminosity_density_lambda[
                self.packet_wvl_range_mask
            ]
            / self.lum_to_flux
        )

    def _calculate_emission_luminosities(self, packets_mode, packet_wvl_range):
        """
        Calculate luminosities for the emission part of SDEC plot.

        Parameters
        ----------
        packets_mode : {'virtual', 'real'}
            Mode of packets to be considered, either real or virtual
        packet_wvl_range : astropy.Quantity
            Wavelength range to restrict the analysis of escaped packets. It
            should be a quantity having units of Angstrom, containing two
            values - lower lambda and upper lambda i.e.
            [lower_lambda, upper_lambda] * u.AA

        Returns
        -------
        luminosities_df : pd.DataFrame
            Dataframe containing luminosities contributed by no interaction,
            only e-scattering and emission with each element present
        elements_present: np.array
            Atomic numbers of the elements with which packets of specified
            wavelength range interacted
        """
        # Calculate masks to be applied on packets data based on packet_wvl_range
        if packet_wvl_range is None:
            self.packet_nu_range_mask = np.ones(
                self.data[packets_mode].packets_df.shape[0], dtype=bool
            )
            self.packet_nu_line_range_mask = np.ones(
                self.data[packets_mode].packets_df_line_interaction.shape[0],
                dtype=bool,
            )
        else:
            packet_nu_range = packet_wvl_range.to("Hz", u.spectral())
            self.packet_nu_range_mask = (
                self.data[packets_mode].packets_df["nus"] < packet_nu_range[0]
            ) & (self.data[packets_mode].packets_df["nus"] > packet_nu_range[1])
            self.packet_nu_line_range_mask = (
                self.data[packets_mode].packets_df_line_interaction["nus"]
                < packet_nu_range[0]
            ) & (
                self.data[packets_mode].packets_df_line_interaction["nus"]
                > packet_nu_range[1]
            )

        # Histogram weights are packet luminosities or flux
        weights = (
            self.data[packets_mode].packets_df["energies"][
                self.packet_nu_range_mask
            ]
            / self.lum_to_flux
        ) / self.data[packets_mode].time_of_simulation

        luminosities_df = pd.DataFrame(index=self.plot_wavelength)

        # Contribution of packets which experienced no interaction ------------
        # Mask to select packets with no interaction
        mask_noint = (
            self.data[packets_mode].packets_df["last_interaction_type"][
                self.packet_nu_range_mask
            ]
            == -1
        )

        # Calculate weighted histogram of packet frequencies for
        # plottable range of frequency bins
        hist_noint = np.histogram(
            self.data[packets_mode].packets_df["nus"][
                self.packet_nu_range_mask
            ][mask_noint],
            bins=self.plot_frequency_bins,
            weights=weights[mask_noint],
            density=False,
        )

        # Convert histogram (luminosity values) to luminosity density lambda
        L_nu_noint = (
            hist_noint[0]
            * u.erg
            / u.s
            / self.data[packets_mode].spectrum_delta_frequency
        )
        L_lambda_noint = L_nu_noint * self.plot_frequency / self.plot_wavelength

        # Save it in df
        luminosities_df["noint"] = L_lambda_noint.value

        # Contribution of packets which only experienced electron scattering ---
        mask_escatter = (
            self.data[packets_mode].packets_df["last_interaction_type"][
                self.packet_nu_range_mask
            ]
            == 1
        ) & (
            self.data[packets_mode].packets_df["last_line_interaction_in_id"][
                self.packet_nu_range_mask
            ]
            == -1
        )
        hist_escatter = np.histogram(
            self.data[packets_mode].packets_df["nus"][
                self.packet_nu_range_mask
            ][mask_escatter],
            bins=self.plot_frequency_bins,
            weights=weights[mask_escatter],
            density=False,
        )

        L_nu_escatter = (
            hist_escatter[0]
            * u.erg
            / u.s
            / self.data[packets_mode].spectrum_delta_frequency
        )
        L_lambda_escatter = (
            L_nu_escatter * self.plot_frequency / self.plot_wavelength
        )
        luminosities_df["escatter"] = L_lambda_escatter.value

        # Group packets_df by atomic number of elements with which packets
        # had their last emission (interaction out)
        packets_df_grouped = (
            self.data[packets_mode]
            .packets_df_line_interaction.loc[self.packet_nu_line_range_mask]
            .groupby(by="last_line_interaction_out_atom")
        )

        # Contribution of each element with which packets interacted ----------
        for atomic_number, group in packets_df_grouped:
            # Histogram of specific element
            hist_el = np.histogram(
                group["nus"],
                bins=self.plot_frequency_bins,
                weights=group["energies"]
                / self.lum_to_flux
                / self.data[packets_mode].time_of_simulation,
            )

            # Convert to luminosity density lambda
            L_nu_el = (
                hist_el[0]
                * u.erg
                / u.s
                / self.data[packets_mode].spectrum_delta_frequency
            )
            L_lambda_el = L_nu_el * self.plot_frequency / self.plot_wavelength

            luminosities_df[atomic_number] = L_lambda_el.value

        # Create an array of the elements with which packets interacted
        elements_present = np.array(list(packets_df_grouped.groups.keys()))

        return luminosities_df, elements_present

    def _calculate_absorption_luminosities(
        self, packets_mode, packet_wvl_range
    ):
        """
        Calculate luminosities for the absorption part of SDEC plot.

        Parameters
        ----------
        packets_mode : {'virtual', 'real'}
            Mode of packets to be considered, either real or virtual
        packet_wvl_range : astropy.Quantity
            Wavelength range to restrict the analysis of escaped packets. It
            should be a quantity having units of Angstrom, containing two
            values - lower lambda and upper lambda i.e.
            [lower_lambda, upper_lambda] * u.AA

        Returns
        -------
        pd.DataFrame
            Dataframe containing luminosities contributed by absorption with
            each element present
        """
        # Calculate masks to be applied on packets data based on packet_wvl_range
        if packet_wvl_range is None:
            self.packet_nu_line_range_mask = np.ones(
                self.data[packets_mode].packets_df_line_interaction.shape[0],
                dtype=bool,
            )
        else:
            packet_nu_range = packet_wvl_range.to("Hz", u.spectral())
            self.packet_nu_line_range_mask = (
                self.data[packets_mode].packets_df_line_interaction[
                    "last_line_interaction_in_nu"
                ]
                < packet_nu_range[0]
            ) & (
                self.data[packets_mode].packets_df_line_interaction[
                    "last_line_interaction_in_nu"
                ]
                > packet_nu_range[1]
            )

        luminosities_df = pd.DataFrame(index=self.plot_wavelength)

        # Group packets_df by atomic number of elements with which packets
        # had their last absorption (interaction in)
        packets_df_grouped = (
            self.data[packets_mode]
            .packets_df_line_interaction.loc[self.packet_nu_line_range_mask]
            .groupby(by="last_line_interaction_in_atom")
        )
        for atomic_number, group in packets_df_grouped:
            # Histogram of specific element
            hist_el = np.histogram(
                group["last_line_interaction_in_nu"],
                bins=self.plot_frequency_bins,
                weights=group["energies"]
                / self.lum_to_flux
                / self.data[packets_mode].time_of_simulation,
            )

            # Convert to luminosity density lambda
            L_nu_el = (
                hist_el[0]
                * u.erg
                / u.s
                / self.data[packets_mode].spectrum_delta_frequency
            )
            L_lambda_el = L_nu_el * self.plot_frequency / self.plot_wavelength

            luminosities_df[atomic_number] = L_lambda_el.value

        return luminosities_df

    def _calculate_photosphere_luminosity(self, packets_mode):
        """
        Calculate blackbody luminosity of the photosphere.

        Parameters
        ----------
        packets_mode : {'virtual', 'real'}
            Mode of packets to be considered, either real or virtual

        Returns
        -------
        astropy.Quantity
            Luminosity density lambda (or Flux) of photosphere (inner boundary
            of TARDIS simulation)
        """
        L_lambda_ph = (
            abb.blackbody_lambda(
                self.plot_wavelength,
                self.data[packets_mode].t_inner,
            )
            * 4
            * np.pi ** 2
            * self.data[packets_mode].r_inner[0] ** 2
            * u.sr
        ).to("erg / (AA s)")

        return L_lambda_ph / self.lum_to_flux

    def generate_plot_mpl(
        self,
        packets_mode="virtual",
        packet_wvl_range=None,
        distance=None,
        show_modeled_spectrum=True,
        ax=None,
        figsize=(10, 7),
        cmapname="jet",
    ):
        """
        Generate Spectral element DEComposition (SDEC) Plot using matplotlib.

        Parameters
        ----------
        packets_mode : {'virtual', 'real'}, optional
            Mode of packets to be considered, either real or virtual. Default
            value is 'virtual'
        packet_wvl_range : astropy.Quantity or None, optional
            Wavelength range to restrict the analysis of escaped packets. It
            should be a quantity having units of Angstrom, containing two
            values - lower lambda and upper lambda i.e.
            [lower_lambda, upper_lambda] * u.AA. Default value is None
        distance : astropy.Quantity or None, optional
            Distance used to calculate flux instead of luminosities in the plot.
            It should have a length unit like m, Mpc, etc. Default value is None
        show_modeled_spectrum : bool, optional
            Whether to show modeled spectrum in SDEC Plot. Default value is
            True
        ax : matplotlib.axes._subplots.AxesSubplot or None, optional
            Axis on which to create plot. Default value is None which will
            create plot on a new figure's axis.
        figsize : tuple, optional
            Size of the matplotlib figure to display. Default value is (10, 7)
        cmapname : str, optional
            Name of matplotlib colormap to be used for showing elements.
            Default value is "jet"

        Returns
        -------
        matplotlib.axes._subplots.AxesSubplot
            Axis on which SDEC Plot is created
        """
        # Calculate data attributes required for plotting
        # and save them in instance itself
        self._calculate_plotting_data(
            packets_mode=packets_mode,
            packet_wvl_range=packet_wvl_range,
            distance=distance,
        )

        if ax is None:
            self.ax = plt.figure(figsize=figsize).add_subplot(111)
        else:
            self.ax = ax

        # Set colormap to be used in elements of emission and absorption plots
        self.cmap = cm.get_cmap(cmapname, self.elements.size)

        self._plot_emission_mpl()
        self._plot_absorption_mpl()

        # Plot modeled spectrum
        if show_modeled_spectrum:
            self.ax.plot(
                self.plot_wavelength,
                self.modeled_spectrum_luminosity,
                "--b",
                label=f"{packets_mode.capitalize()} Spectrum",
                linewidth=1,
            )

        # Plot photosphere
        self.ax.plot(
            self.plot_wavelength,
            self.photosphere_luminosity,
            "--r",
            label="Blackbody Photosphere",
        )

        self._show_colorbar_mpl()

        # Set legends and labels
        self.ax.legend(fontsize=12)
        self.ax.set_xlabel(r"Wavelength $[\AA]$", fontsize=12)
        if distance:  # Set y-axis label for flux
            self.ax.set_ylabel(
                r"$F_{\lambda}$ [erg/s/cm^{2}/$\AA$]", fontsize=12
            )
        else:  # Set y-axis label for luminosity
            self.ax.set_ylabel(r"$L_{\lambda}$ [erg/s/$\AA$]", fontsize=12)

        return plt.gca()

    def _plot_emission_mpl(self):
        """Plot emission part of the SDEC Plot using matplotlib."""
        # To create stacked area chart in matplotlib, we will start with zero
        # lower level and will keep adding luminosities to it (upper level)
        lower_level = np.zeros(self.emission_luminosities_df.shape[0])
        upper_level = (
            lower_level + self.emission_luminosities_df.noint.to_numpy()
        )

        self.ax.fill_between(
            self.plot_wavelength,
            lower_level,
            upper_level,
            color="black",
            label="No interaction",
        )

        lower_level = upper_level
        upper_level = (
            lower_level + self.emission_luminosities_df.escatter.to_numpy()
        )

        self.ax.fill_between(
            self.plot_wavelength,
            lower_level,
            upper_level,
            color="grey",
            label="Electron Scatter Only",
        )

        elements_z = self.emission_luminosities_df.columns[2:].to_list()
        nelements = len(elements_z)

        # Contribution from each element
        for i, atomic_number in enumerate(elements_z):
            lower_level = upper_level
            upper_level = (
                lower_level
                + self.emission_luminosities_df[atomic_number].to_numpy()
            )

            self.ax.fill_between(
                self.plot_wavelength,
                lower_level,
                upper_level,
                color=self.cmap(i / nelements),
                cmap=self.cmap,
                linewidth=0,
            )

    def _plot_absorption_mpl(self):
        """Plot absorption part of the SDEC Plot using matplotlib."""
        lower_level = np.zeros(self.absorption_luminosities_df.shape[0])

        elements_z = self.absorption_luminosities_df.columns.to_list()
        for i, atomic_number in enumerate(elements_z):
            # To plot absorption part along -ve X-axis, we will start with
            # zero upper level and keep subtracting luminosities to it (lower
            # level) - fill from upper to lower level
            upper_level = lower_level
            lower_level = (
                upper_level
                - self.absorption_luminosities_df[atomic_number].to_numpy()
            )

            self.ax.fill_between(
                self.plot_wavelength,
                upper_level,
                lower_level,
                color=self.cmap(i / len(elements_z)),
                cmap=self.cmap,
                linewidth=0,
            )

    def _show_colorbar_mpl(self):
        """Show matplotlib colorbar with labels of elements mapped to colors."""
        color_values = [
            self.cmap(i / self.elements.size) for i in range(self.elements.size)
        ]
        custcmap = clr.ListedColormap(color_values)
        norm = clr.Normalize(vmin=0, vmax=self.elements.size)
        mappable = cm.ScalarMappable(norm=norm, cmap=custcmap)
        mappable.set_array(np.linspace(1, self.elements.size + 1, 256))
        cbar = plt.colorbar(mappable, ax=self.ax)

        bounds = np.arange(self.elements.size) + 0.5
        cbar.set_ticks(bounds)

        elements_name = [
            atomic_number2element_symbol(atomic_num)
            for atomic_num in self.elements
        ]
        cbar.set_ticklabels(elements_name)

    def generate_plot_ply(
        self,
        packets_mode="virtual",
        packet_wvl_range=None,
        distance=None,
        observed_spectrum=None,
        show_modeled_spectrum=True,
        fig=None,
        graph_height=600,
        cmapname="jet",
    ):
        """
        Generate interactive Spectral element DEComposition (SDEC) Plot using plotly.

        Parameters
        ----------
        packets_mode : {'virtual', 'real'}, optional
            Mode of packets to be considered, either real or virtual. Default
            value is 'virtual'
        packet_wvl_range : astropy.Quantity or None, optional
            Wavelength range to restrict the analysis of escaped packets. It
            should be a quantity having units of Angstrom, containing two
            values - lower lambda and upper lambda i.e.
            [lower_lambda, upper_lambda] * u.AA. Default value is None
        distance : astropy.Quantity or None, optional
            Distance used to calculate flux instead of luminosities in the plot.
            It should have a length unit like m, Mpc, etc. Default value is None
        show_modeled_spectrum : bool, optional
            Whether to show modeled spectrum in SDEC Plot. Default value is
            True
        fig : plotly.graph_objs._figure.Figure or None, optional
            Figure object on which to create plot. Default value is None which
            will create plot on a new Figure object.
        graph_height : int, optional
            Height (in px) of the plotly graph to display. Default value is 600
        cmapname : str, optional
            Name of the colormap to be used for showing elements.
            Default value is "jet"

        Returns
        -------
        plotly.graph_objs._figure.Figure
            Figure object on which SDEC Plot is created
        """
        # Calculate data attributes required for plotting
        # and save them in instance itself
        self._calculate_plotting_data(
            packets_mode=packets_mode,
            packet_wvl_range=packet_wvl_range,
            distance=distance,
        )

        if fig is None:
            self.fig = go.Figure()
        else:
            self.fig = fig

        # Set colormap to be used in elements of emission and absorption plots
        self.cmap = cm.get_cmap(cmapname, self.elements.size)

        self._plot_emission_ply()
        self._plot_absorption_ply()

        # Plot modeled spectrum
        if show_modeled_spectrum:
            self.fig.add_trace(
                go.Scatter(
                    x=self.plot_wavelength,
                    y=self.modeled_spectrum_luminosity,
                    mode="lines",
                    line=dict(
                        color="blue",
                        width=1,
                    ),
                    name=f"{packets_mode.capitalize()} Spectrum",
                )
            )

        # Plot photosphere
        self.fig.add_trace(
            go.Scatter(
                x=self.plot_wavelength,
                y=self.photosphere_luminosity,
                mode="lines",
                line=dict(width=1.5, color="red", dash="dash"),
                name="Blackbody Photosphere",
            )
        )

        self._show_colorbar_ply()

        # Set label and other layout options
        xlabel = pu.axis_label_in_latex("Wavelength", u.AA)
        if distance:  # Set y-axis label for flux
            ylabel = pu.axis_label_in_latex(
                "F_{\\lambda}", u.Unit("erg/(s cm**2 AA)"), only_text=False
            )
        else:  # Set y-axis label for luminosity
            ylabel = pu.axis_label_in_latex(
                "L_{\\lambda}", u.Unit("erg/(s AA)"), only_text=False
            )
        self.fig.update_layout(
            xaxis=dict(
                title=xlabel,
                exponentformat="none",
            ),
            yaxis=dict(title=ylabel, exponentformat="e"),
            height=graph_height,
        )

        return self.fig

    @staticmethod
    def to_rgb255_string(color_tuple):
        """
        Convert a matplotlib RGBA tuple to a generic RGB 255 string.

        Parameters
        ----------
        color_tuple : tuple
            Matplotlib RGBA tuple of float values in closed interval [0, 1]

        Returns
        -------
        str
            RGB string of format rgb(r,g,b) where r,g,b are integers between
            0 and 255 (both inclusive)
        """
        color_tuple_255 = tuple([int(x * 255) for x in color_tuple[:3]])
        return f"rgb{color_tuple_255}"

    def _plot_emission_ply(self):
        """Plot emission part of the SDEC Plot using plotly."""
        # By specifying a common stackgroup, plotly will itself add up
        # luminosities, in order, to created stacked area chart
        self.fig.add_trace(
            go.Scatter(
                x=self.emission_luminosities_df.index,
                y=self.emission_luminosities_df.noint,
                mode="none",
                name="No interaction",
                fillcolor="black",
                stackgroup="emission",
            )
        )

        self.fig.add_trace(
            go.Scatter(
                x=self.emission_luminosities_df.index,
                y=self.emission_luminosities_df.escatter,
                mode="none",
                name="Electron Scatter Only",
                fillcolor="grey",
                stackgroup="emission",
            )
        )

        elements_z = self.emission_luminosities_df.columns[2:]
        nelements = len(elements_z)

        for i, atomic_num in enumerate(elements_z):
            self.fig.add_trace(
                go.Scatter(
                    x=self.emission_luminosities_df.index,
                    y=self.emission_luminosities_df[atomic_num],
                    mode="none",
                    name=atomic_number2element_symbol(atomic_num),
                    fillcolor=self.to_rgb255_string(self.cmap(i / nelements)),
                    stackgroup="emission",
                    showlegend=False,
                )
            )

    def _plot_absorption_ply(self):
        """Plot absorption part of the SDEC Plot using plotly."""
        elements_z = self.absorption_luminosities_df.columns
        nelements = len(elements_z)

        for i, atomic_num in enumerate(elements_z):
            self.fig.add_trace(
                go.Scatter(
                    x=self.absorption_luminosities_df.index,
                    # to plot absorption luminosities along negative y-axis
                    y=self.absorption_luminosities_df[atomic_num] * -1,
                    mode="none",
                    name=atomic_number2element_symbol(atomic_num),
                    fillcolor=self.to_rgb255_string(self.cmap(i / nelements)),
                    stackgroup="absorption",
                    showlegend=False,
                )
            )

    def _show_colorbar_ply(self):
        """Show plotly colorbar with labels of elements mapped to colors."""
        # Interpolate [0, 1] range to create bins equal to number of elements
        colorscale_bins = np.linspace(0, 1, num=self.elements.size + 1)

        # Create a categorical colorscale [a list of (reference point, color)]
        # by mapping same reference points (excluding 1st and last bin edge)
        # twice in a row (https://plotly.com/python/colorscales/#constructing-a-discrete-or-discontinuous-color-scale)
        categorical_colorscale = []
        for i in range(self.elements.size):
            color = self.to_rgb255_string(self.cmap(colorscale_bins[i]))
            categorical_colorscale.append((colorscale_bins[i], color))
            categorical_colorscale.append((colorscale_bins[i + 1], color))

        coloraxis_options = dict(
            colorscale=categorical_colorscale,
            showscale=True,
            cmin=0,
            cmax=self.elements.size,
            colorbar=dict(
                title="Elements",
                tickvals=np.arange(0, self.elements.size) + 0.5,
                ticktext=[
                    atomic_number2element_symbol(atomic_num)
                    for atomic_num in self.elements
                ],
                # to change length and position of colorbar
                len=0.75,
                yanchor="top",
                y=0.75,
            ),
        )

        # Plot an invisible one point scatter trace, to make colorbar show up
        scatter_point_idx = pu.get_mid_point_idx(self.plot_wavelength)
        self.fig.add_trace(
            go.Scatter(
                x=self.plot_wavelength[scatter_point_idx],
                y=[0],
                mode="markers",
                showlegend=False,
                hoverinfo="skip",
                marker=dict(color=[0], opacity=0, **coloraxis_options),
            )
        )
