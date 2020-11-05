"""
Interface to generate Kromer Plot for TARDIS simulation models.

Kromer Plot is a spectral diagnostics plot similar to those originally
proposed by M. Kromer (see, for example, Kromer et al. 2013, figure 4).
"""
import numpy as np
import pandas as pd
import pyne
import astropy.units as u
import astropy.modeling.blackbody as abb

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm


class KromerData:
    """The data of simulation model which is used by Kromer Plot."""

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
        Initialize the KromerData with required properties of simulation model.

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
        )
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
        Create an instance of KromerData from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation
        packets_mode : {'virtual', 'real'}
            Mode of packets to be considered, either real or virtual

        Returns
        -------
        KromerData
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
        Create an instance of KromerData from a simulation HDF file.

        Parameters
        ----------
        hdf_fpath : str
            Valid path to the HDF file where simulation is saved
        packets_mode : {'virtual', 'real'}
            Mode of packets to be considered, either real or virtual

        Returns
        -------
        KromerData
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


class KromerPlotter:
    """Plotting interface to generate Kromer Plot for a simulation model."""

    def __init__(self, data):
        """
        Initialize the KromerPlotter with required data of simulation model.

        Parameters
        ----------
        data : dict of KromerData
            Dictionary to store data required for Kromer plot, for both packet
            modes i.e. real and virtual
        """
        self.data = data

    @classmethod
    def from_simulation(cls, sim):
        """
        Create an instance of KromerPlotter from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS Simulation object produced by running a simulation

        Returns
        -------
        KromerPlotter
        """
        return cls(
            dict(
                virtual=KromerData.from_simulation(sim, "virtual"),
                real=KromerData.from_simulation(sim, "real"),
            )
        )

    @classmethod
    def from_hdf(cls, hdf_fpath):
        """
        Create an instance of KromerPlotter from a simulation HDF file.

        Parameters
        ----------
        hdf_fpath : str
            Valid path to the HDF file where simulation is saved

        Returns
        -------
        KromerPlotter
        """
        return cls(
            dict(
                virtual=KromerData.from_hdf(hdf_fpath, "virtual"),
                real=KromerData.from_hdf(hdf_fpath, "real"),
            )
        )

    def generate_plot(
        self,
        packets_mode="virtual",
        packet_wvl_range=None,
        distance=None,
        observed_spectrum=None,
        show_modeled_spectrum=True,
        ax=None,
        figsize=(10, 7),
        cmapname="jet",
    ):
        # TODO: Read observed_spectrum in ascii or csv - 2 column format

        if packets_mode not in ["virtual", "real"]:
            raise ValueError(
                "Invalid value passed to packets_mode. Only "
                "allowed values are 'virtual' or 'real'"
            )

        if distance is None:
            self.lum_to_flux = 1  # so that this factor will have no effect
        else:
            self.lum_to_flux = 4.0 * np.pi * (distance.to("cm")) ** 2

        # Bin edges
        bins = self.data[packets_mode].spectrum_frequency_bins

        # Wavelengths
        wvl = self.data[packets_mode].spectrum_wavelength

        emission_luminosities_df = self._calculate_emission_luminosities(
            packets_mode=packets_mode,
            bins=bins,
            wvl=wvl,
            packet_wvl_range=packet_wvl_range,
        )
        absorption_luminosities_df = self._calculate_absorption_luminosities(
            packets_mode=packets_mode,
            bins=bins,
            wvl=wvl,
            packet_wvl_range=packet_wvl_range,
        )
        photosphere_luminosity = self._calculate_photosphere_luminosity(
            packets_mode=packets_mode
        )

        # Plotting ------------------------------------------------
        if ax is None:
            self.ax = plt.figure(figsize=figsize).add_subplot(111)
        else:
            self.ax = ax

        # Plot virtual spectrum
        if show_modeled_spectrum:
            self.ax.plot(
                self.data[packets_mode].spectrum_wavelength,
                self.data[packets_mode].spectrum_luminosity_density_lambda
                / self.lum_to_flux,
                "--b",
                label=f"{packets_mode.capitalize()} Spectrum",
                # ds="steps-pre", # no need to make it look histogram
                linewidth=1,
            )

        self._plot_emission(
            luminosities_df=emission_luminosities_df, cmapname=cmapname
        )
        self._plot_absorption(absorption_luminosities_df)

        # Plot Photosphere luminosity
        self.ax.plot(
            self.data[packets_mode].spectrum_wavelength,
            photosphere_luminosity / self.lum_to_flux,
            "--r",
            label="Blackbody Photosphere",
        )

        self.ax.legend(fontsize=12)
        self.ax.set_xlabel(r"Wavelength $(\AA)$", fontsize=15)
        if distance:  # Set y-axis label for flux
            self.ax.set_ylabel(
                r"$F_{\lambda}$ (erg/s/$cm^{2}/\AA$)", fontsize=15
            )
        else:  # Set y-axis label for luminosity
            self.ax.set_ylabel(r"$L$ (erg/s/$\AA$)", fontsize=15)

        return plt.gca()

    def _calculate_emission_luminosities(
        self,
        packets_mode,
        bins,
        wvl,
        packet_wvl_range,
    ):
        if packet_wvl_range is None:
            packet_nu_range_mask = np.ones(
                self.data[packets_mode].packets_df.shape[0], dtype=bool
            )
            packet_nu_line_range_mask = np.ones(
                self.data[packets_mode].packets_df_line_interaction.shape[0],
                dtype=bool,
            )
        else:
            packet_nu_range = packet_wvl_range.to("Hz", u.spectral())
            # print(packet_nu_range)
            packet_nu_range_mask = (
                self.data[packets_mode].packets_df["nus"] < packet_nu_range[0]
            ) & (self.data[packets_mode].packets_df["nus"] > packet_nu_range[1])
            packet_nu_line_range_mask = (
                self.data[packets_mode].packets_df_line_interaction["nus"]
                < packet_nu_range[0]
            ) & (
                self.data[packets_mode].packets_df_line_interaction["nus"]
                > packet_nu_range[1]
            )

        # weights are packet luminosities or flux
        weights = (
            self.data[packets_mode].packets_df["energies"][packet_nu_range_mask]
            / self.lum_to_flux
        ) / self.data[packets_mode].time_of_simulation

        luminosities_df = pd.DataFrame(index=wvl)

        # No interaction contribution
        # mask_noint selects packets with no interaction
        mask_noint = (
            self.data[packets_mode].packets_df["last_interaction_type"][
                packet_nu_range_mask
            ]
            == -1
        )
        hist_noint = np.histogram(
            self.data[packets_mode].packets_df["nus"][packet_nu_range_mask][
                mask_noint
            ],
            bins=bins,
            weights=weights[mask_noint],
            density=False,
        )

        # We convert histogram values to luminosity density lambda.
        L_nu_noint = (
            hist_noint[0]
            * u.erg
            / u.s
            / self.data[packets_mode].spectrum_delta_frequency
        )
        L_lambda_noint = (
            L_nu_noint
            * self.data[packets_mode].spectrum_frequency
            / self.data[packets_mode].spectrum_wavelength
        )
        # Save it in df
        luminosities_df["noint"] = L_lambda_noint.value

        # Electron scattering contribution
        # mask_escatter selects packets that ONLY experience
        # electron scattering.
        mask_escatter = (
            self.data[packets_mode].packets_df["last_interaction_type"][
                packet_nu_range_mask
            ]
            == 1
        ) & (
            self.data[packets_mode].packets_df["last_line_interaction_in_id"][
                packet_nu_range_mask
            ]
            == -1
        )
        hist_escatter = np.histogram(
            self.data[packets_mode].packets_df["nus"][packet_nu_range_mask][
                mask_escatter
            ],
            bins=bins,
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
            L_nu_escatter
            * self.data[packets_mode].spectrum_frequency
            / self.data[packets_mode].spectrum_wavelength
        )
        luminosities_df["escatter"] = L_lambda_escatter.value

        # Groupby packet dataframe by last atom interaction out
        g = (
            self.data[packets_mode]
            .packets_df_line_interaction.loc[packet_nu_line_range_mask]
            .groupby(by="last_line_interaction_out_atom")
        )

        # Contribution from each element
        for atomic_number, group in g:
            # histogram of specific element
            hist_el = np.histogram(
                group["nus"],
                bins=bins,
                weights=group["energies"]
                / self.data[packets_mode].time_of_simulation,
            )

            # Convert to luminosity density lambda
            L_nu_el = (
                hist_el[0]
                * u.erg
                / u.s
                / self.data[packets_mode].spectrum_delta_frequency
            )
            L_lambda_el = (
                L_nu_el
                * self.data[packets_mode].spectrum_frequency
                / self.data[packets_mode].spectrum_wavelength
            )

            luminosities_df[atomic_number] = L_lambda_el.value

        return luminosities_df

    def _plot_emission(self, luminosities_df, cmapname):
        wavelength = luminosities_df.index.to_numpy()

        lower_level = np.zeros(luminosities_df.shape[0])
        upper_level = lower_level + luminosities_df.noint.to_numpy()

        self.ax.fill_between(
            wavelength,
            lower_level,
            upper_level,
            # step="pre",
            color="k",
            label="No interaction",
        )

        lower_level = upper_level
        upper_level = lower_level + luminosities_df.escatter.to_numpy()

        self.ax.fill_between(
            wavelength,
            lower_level,
            upper_level,
            # step="pre",
            color="grey",
            label="Electron Scatter Only",
        )

        # Set up color map
        elements_z = luminosities_df.columns[2:].to_list()
        nelements = len(elements_z)

        # Save color map for later use
        self.cmap = cm.get_cmap(cmapname, nelements)

        # Contribution from each element
        for i, atomic_number in enumerate(elements_z):
            lower_level = upper_level
            upper_level = (
                lower_level + luminosities_df[atomic_number].to_numpy()
            )

            self.ax.fill_between(
                wavelength,
                lower_level,
                upper_level,
                # step="pre",
                color=self.cmap(i / nelements),
                cmap=self.cmap,
                linewidth=0,
            )

        values = [self.cmap(i / nelements) for i in range(nelements)]
        custcmap = matplotlib.colors.ListedColormap(values)
        bounds = np.arange(nelements) + 0.5
        norm = matplotlib.colors.Normalize(vmin=0, vmax=nelements)
        mappable = cm.ScalarMappable(norm=norm, cmap=custcmap)
        mappable.set_array(np.linspace(1, nelements + 1, 256))

        elements_name = [pyne.nucname.name(el) for el in elements_z]
        labels = elements_name

        cbar = plt.colorbar(mappable, ax=self.ax)
        cbar.set_ticks(bounds)
        cbar.set_ticklabels(labels)

    def _calculate_absorption_luminosities(
        self, packets_mode, bins, wvl, packet_wvl_range
    ):

        # Set up a wavelength range to only analyze emitted packets in that range.
        # packet_wvl_range should be a list [lower_lambda, upper_lambda]*u.angstrom
        # Notice that we convert to frequency, where the lower_lambda is a larger
        # frequency.
        if packet_wvl_range is None:
            packet_nu_line_range_mask = np.ones(
                self.data[packets_mode].packets_df_line_interaction.shape[0],
                dtype=bool,
            )
        else:
            packet_nu_range = packet_wvl_range.to("Hz", u.spectral())
            packet_nu_line_range_mask = (
                self.data[packets_mode].packets_df_line_interaction["nus"]
                < packet_nu_range[0]
            ) & (
                self.data[packets_mode].packets_df_line_interaction["nus"]
                > packet_nu_range[1]
            )

        luminosities_df = pd.DataFrame(index=wvl)

        # Groupby packet dataframe by last atom interaction in
        g_abs = (
            self.data[packets_mode]
            .packets_df_line_interaction.loc[packet_nu_line_range_mask]
            .groupby(by="last_line_interaction_in_atom")
        )
        for atomic_number, group in g_abs:
            hist_el = np.histogram(
                group["last_line_interaction_in_nu"],
                bins=bins,
                weights=group["energies"]
                / self.data[packets_mode].time_of_simulation,
            )
            L_nu_el = (
                hist_el[0]
                * u.erg
                / u.s
                / self.data[packets_mode].spectrum_delta_frequency
            )
            L_lambda_el = (
                L_nu_el
                * self.data[packets_mode].spectrum_frequency
                / self.data[packets_mode].spectrum_wavelength
            )

            luminosities_df[atomic_number] = L_lambda_el.value

        return luminosities_df

    def _plot_absorption(self, luminosities_df):
        wavelength = luminosities_df.index.to_numpy()
        lower_level = np.zeros(luminosities_df.shape[0])

        elements_z = luminosities_df.columns.to_list()
        for i, atomic_number in enumerate(elements_z):
            # Fill from upper to lower level, moving along -ve x-axis
            upper_level = lower_level
            lower_level = (
                upper_level - luminosities_df[atomic_number].to_numpy()
            )

            self.ax.fill_between(
                wavelength,
                upper_level,
                lower_level,
                # step="pre",
                color=self.cmap(i / len(elements_z)),
                cmap=self.cmap,
                linewidth=0,
            )

    def _calculate_photosphere_luminosity(self, packets_mode):
        """
        Plots the blackbody luminosity density from the inner
        boundary of the TARDIS simulation.

        """
        return (
            abb.blackbody_lambda(
                self.data[packets_mode].spectrum_wavelength,
                self.data[packets_mode].t_inner,
            )
            * 4
            * np.pi ** 2
            * self.data[packets_mode].r_inner[0] ** 2
            * u.sr
        ).to("erg / (AA s)")
