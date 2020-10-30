"""A simple plotting tool to create spectral diagnostics plots similar to those
originally proposed by M. Kromer (see, for example, Kromer et al. 2013, figure
4).
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
        spectrum_frequency,  # stores _frequency not frequency
        spectrum_luminosity_density_lambda,
        spectrum_wavelength,
        t_inner,
        time_of_simulation,
    ):
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
        self.spectrum_frequency = spectrum_frequency
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
                spectrum_frequency=sim.runner.spectrum_virtual._frequency,
                spectrum_luminosity_density_lambda=sim.runner.spectrum_virtual.luminosity_density_lambda,
                spectrum_wavelength=sim.runner.spectrum_virtual.wavelength,
                t_inner=t_inner,
                time_of_simulation=time_of_simulation,
            )

        elif packets_mode == "real":
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
                spectrum_frequency=sim.runner.spectrum._frequency,
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
        with pd.HDFStore(hdf_fpath, "r") as hdf:
            lines_df = (
                hdf["/simulation/plasma/lines"]
                .reset_index()
                .set_index("line_id")
            )
            r_inner = u.Quantity(
                hdf["/simulation/model/r_inner"].to_numpy(), "cm"
            )  # Convert series to array to construct quantity from it
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
                    spectrum_frequency=u.Quantity(
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
                    spectrum_frequency=u.Quantity(
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
    def __init__(self, data):
        self.data = data

    @classmethod
    def from_simulation(cls, sim):
        return cls(
            dict(
                virtual=KromerData.from_simulation(sim, "virtual"),
                real=KromerData.from_simulation(sim, "real"),
            )
        )

    @classmethod
    def from_hdf(cls, hdf_fpath):
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

        # Read ascii or csv - 2 column format
        if packets_mode not in ["virtual", "real"]:
            raise ValueError(
                "Invalid value passed to packets_mode. Only "
                "allowed values are 'virtual' or 'real'"
            )

        if distance is None:
            self.lum_to_flux = 1  # so that this factor will have no effect
        else:
            self.lum_to_flux = 4.0 * np.pi * (distance.to("cm")) ** 2

        if ax is None:
            self.ax = plt.figure(figsize=figsize).add_subplot(111)
        else:
            self.ax = ax

        # Bin edges
        bins = self.data[packets_mode].spectrum_frequency

        # Wavelengths
        wvl = self.data[packets_mode].spectrum_wavelength

        self._plot_emission(
            packets_mode=packets_mode,
            bins=bins,
            wvl=wvl,
            packet_wvl_range=packet_wvl_range,
            cmapname=cmapname,
            show_modeled_spectrum=show_modeled_spectrum,
        )
        self._plot_absorption(
            packets_mode=packets_mode,
            bins=bins,
            wvl=wvl,
            packet_wvl_range=packet_wvl_range,
        )
        self._plot_photosphere(packets_mode=packets_mode)

        self.ax.legend(fontsize=12)
        self.ax.set_xlabel(r"Wavelength $(\AA)$", fontsize=15)
        if distance:  # Set y-axis label for flux
            self.ax.set_ylabel(
                r"$F_{\lambda}$ (erg/s/$cm^{2}/\AA$)", fontsize=15
            )
        else:  # Set y-axis label for luminosity
            self.ax.set_ylabel(r"$L$ (erg/s/$\AA$)", fontsize=15)

        return plt.gca()

    def _plot_emission(
        self,
        packets_mode,
        bins,
        wvl,
        packet_wvl_range,
        cmapname,
        show_modeled_spectrum,
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
        hist = np.histogram(
            self.data[packets_mode].packets_df["nus"][packet_nu_range_mask],
            bins=bins,
            weights=weights,
            density=False,
        )

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

        # No Scattering Contribution
        # We convert histogram values to luminosity density lambda.
        lower_level = np.zeros(len(hist_noint[1][1:]))
        L_nu_noint = (
            hist_noint[0]
            * u.erg
            / u.s
            / self.data[packets_mode].spectrum_delta_frequency
        )
        L_lambda_noint = (
            L_nu_noint
            * self.data[packets_mode].spectrum_frequency[:-1]
            / self.data[packets_mode].spectrum_wavelength
        )

        self.ax.fill_between(
            wvl,
            lower_level,
            L_lambda_noint.value,
            # step="pre",
            color="k",
            label="No interaction",
        )
        lower_level = L_lambda_noint.value

        # Only Electron Scattering
        L_nu_escatter = (
            hist_escatter[0]
            * u.erg
            / u.s
            / self.data[packets_mode].spectrum_delta_frequency
        )
        L_lambda_escatter = (
            L_nu_escatter
            * self.data[packets_mode].spectrum_frequency[:-1]
            / self.data[packets_mode].spectrum_wavelength
        )

        self.ax.fill_between(
            wvl,
            lower_level,
            lower_level + L_lambda_escatter.value,
            # step="pre",
            color="grey",
            label="Electron Scatter Only",
        )
        lower_level = lower_level + L_lambda_escatter.value

        # Groupby packet dataframe by last atom interaction out
        g = (
            self.data[packets_mode]
            .packets_df_line_interaction.loc[packet_nu_line_range_mask]
            .groupby(by="last_line_interaction_out_atom")
        )

        # Set up color map
        elements_z = g.groups.keys()
        elements_name = [pyne.nucname.name(el) for el in elements_z]
        nelements = len(elements_z)

        # Save color map for later use
        self.cmap = cm.get_cmap(cmapname, nelements)

        values = [self.cmap(i / nelements) for i in range(nelements)]
        custcmap = matplotlib.colors.ListedColormap(values)
        bounds = np.arange(nelements) + 0.5
        norm = matplotlib.colors.Normalize(vmin=0, vmax=nelements)
        mappable = cm.ScalarMappable(norm=norm, cmap=custcmap)

        mappable.set_array(np.linspace(1, nelements + 1, 256))
        labels = elements_name

        # Contribution from each element
        for ind, groupkey in enumerate(elements_z):
            # select subgroup of packet dataframe for specific element.
            group = g.get_group(groupkey)

            # histogram specific element.
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
                * self.data[packets_mode].spectrum_frequency[:-1]
                / self.data[packets_mode].spectrum_wavelength
            )

            self.ax.fill_between(
                wvl,
                lower_level,
                lower_level + L_lambda_el.value,
                # step="pre",
                color=self.cmap(ind / len(elements_z)),
                cmap=self.cmap,
                linewidth=0,
            )
            lower_level = lower_level + L_lambda_el.value

        # Colorbar and Legend
        # self.ax.legend(fontsize=20)
        # self.ax.set_xlabel('Wavelength $(\AA)$', fontsize=20)
        # self.ax.set_ylabel('$L$ (erg/s/$\AA$)', fontsize=20)

        cbar = plt.colorbar(mappable, ax=self.ax)
        cbar.set_ticks(bounds)
        cbar.set_ticklabels(labels)

        return

    def _plot_absorption(self, packets_mode, bins, wvl, packet_wvl_range):

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

        abs_lower_level = np.zeros(len(wvl))

        # Groupby packet dataframe by last atom interaction in
        g_abs = (
            self.data[packets_mode]
            .packets_df_line_interaction.loc[packet_nu_line_range_mask]
            .groupby(by="last_line_interaction_in_atom")
        )
        elements_z = g_abs.groups.keys()
        for ind, groupkey in enumerate(elements_z):
            group = g_abs.get_group(groupkey)
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
                * self.data[packets_mode].spectrum_frequency[:-1]
                / self.data[packets_mode].spectrum_wavelength
            )

            self.ax.fill_between(
                wvl,
                abs_lower_level,
                abs_lower_level - L_lambda_el.value,
                # step="pre",
                color=self.cmap(ind / len(elements_z)),
                cmap=self.cmap,
                linewidth=0,
            )
            abs_lower_level = abs_lower_level - L_lambda_el.value

        return

    def _plot_photosphere(self, packets_mode):
        """
        Plots the blackbody luminosity density from the inner
        boundary of the TARDIS simulation.

        """
        Lph = (
            abb.blackbody_lambda(
                self.data[packets_mode].spectrum_wavelength,
                self.data[packets_mode].t_inner,
            )
            * 4
            * np.pi ** 2
            * self.data[packets_mode].r_inner[0] ** 2
            * u.sr
        ).to("erg / (AA s)")
        self.ax.plot(
            self.data[packets_mode].spectrum_wavelength,
            Lph / self.lum_to_flux,
            "--r",
            label="Blackbody Photosphere",
        )
        return
