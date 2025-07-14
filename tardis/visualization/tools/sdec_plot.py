"""
Spectral element DEComposition (SDEC) Plot for TARDIS simulation models.

This plot is a spectral diagnostics plot similar to those originally
proposed by M. Kromer (see, for example, Kromer et al. 2013, figure 4).
"""

import logging

import astropy.units as u
import matplotlib.cm as cm
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from astropy.modeling.models import BlackBody

from tardis.util.base import (
    atomic_number2element_symbol,
    int_to_roman,
)
from tardis.visualization import plot_util as pu

logger = logging.getLogger(__name__)


class SDECPlotter:
    """
    Plotting interface for Spectral element DEComposition (SDEC) Plot.

    It performs necessary calculations to generate SDEC Plot for a simulation
    model, and allows to plot it in matplotlib and plotly.
    """

    def __init__(self):
        """
        Initialize the SDECPlotter with required data of simulation model.
        """
        self.packet_data = {
            "real": {"packets_df": None, "packets_df_line_interaction": None},
            "virtual": {
                "packets_df": None,
                "packets_df_line_interaction": None,
            },
        }
        self.spectrum = {"virtual": None, "real": None}
        self.t_inner = None
        self.r_inner = None
        self.time_of_simulation = None
        self._default_scatter_kwargs = {
            "mode": "none",
            "hovertemplate": "(%{x:.2f}, %{y:.3g})",
        }

        self._predefined_traces = {
            "emission": {
                "noint": {"name": "No interaction", "fillcolor": "#4C4C4C"},
                "escatter": {
                    "name": "Electron Scatter Only",
                    "fillcolor": "#8F8F8F",
                    "hoverlabel": {"namelength": -1},
                },
                "other": {"name": "Other elements", "fillcolor": "#C2C2C2"},
            },
            "absorption": {
                "other": {"name": "Other elements", "fillcolor": "#C2C2C2"},
            },
        }

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
        plotter = cls()
        plotter.t_inner = sim.simulation_state.packet_source.temperature
        plotter.r_inner = sim.simulation_state.geometry.r_inner_active
        plotter.time_of_simulation = (
            sim.transport.transport_state.packet_collection.time_of_simulation
            * u.s
        )

        modes = ["real"]
        if sim.transport.transport_state.virt_logging:
            modes.append("virtual")

        for mode in modes:
            plotter.spectrum[mode] = pu.get_spectrum_data(mode, sim)
            plotter.packet_data[mode] = pu.extract_and_process_packet_data(
                sim, mode
            )

        return plotter

    @classmethod
    def from_hdf(cls, hdf_fpath):
        """
        Create an instance of SDECPlotter from a simulation HDF file.

        Parameters
        ----------
        hdf_fpath : str
            Valid path to the HDF file where simulation is saved
        packets_mode : {'virtual', 'real'}, optional
            Mode of packets to be considered (default: 'virtual')

        Returns
        -------
        SDECPlotter
        """
        plotter = cls()
        with pd.HDFStore(hdf_fpath, "r") as hdf:
            plotter.r_inner = u.Quantity(
                hdf["/simulation/simulation_state/r_inner"].to_numpy(), "cm"
            )
            plotter.t_inner = u.Quantity(
                hdf["/simulation/simulation_state/scalars"].t_inner, "K"
            )
            transport_state_scalars = hdf["/simulation/transport/transport_state/scalars"]
            plotter.time_of_simulation = u.Quantity(
                transport_state_scalars.time_of_simulation,
                "s",
            )

            has_virtual = bool(getattr(transport_state_scalars, "virt_logging", False))
            modes = ["real"] + (["virtual"] if has_virtual else [])

            for mode in modes:
                plotter.spectrum[mode] = pu.extract_spectrum_data_hdf(hdf, mode)
                plotter.packet_data[mode] = pu.extract_and_process_packet_data_hdf(hdf, mode)

        return plotter

    def _parse_species_list(self, species_list):
        """
        Parse user requested species list and create list of species ids to be used.

        Parameters
        ----------
        species_list : list of species to plot
            List of species (e.g. Si II, Ca II, etc.) that the user wants to show as unique colours.
            Species can be given as an ion (e.g. Si II), an element (e.g. Si), a range of ions
            (e.g. Si I - V), or any combination of these (e.g. species_list = [Si II, Fe I-V, Ca])

        """
        if species_list is not None:
            (
                species_mapped_tuples,
                requested_species_ids_tuples,
                keep_colour,
                full_species_list,
            ) = pu.parse_species_list_util(species_list)
            self._full_species_list = full_species_list
            self._species_list = [
                atomic_num * 100 + ion_num
                for atomic_num, ion_num in requested_species_ids_tuples
            ]

            self._species_mapped = {
                (k[0] * 100 + k[1]): [v[0] * 100 + v[1] for v in values]
                for k, values in species_mapped_tuples.items()
            }
            self._keep_colour = keep_colour
        else:
            self._full_species_list = None
            self._species_list = None
            self._species_mapped = None
            self._keep_colour = None

    def _calculate_plotting_data(
        self, packets_mode, packet_wvl_range, distance, nelements
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
            Distance used to calculate flux instead of luminosity in the plot.
            It should have a length unit like m, Mpc, etc.
        nelements: int
            Number of elements to include in plot. Determined by the
            largest contribution to the total luminosity absorbed and emitted.
            Other elements are shown in silver. Default value is
            None, which displays all elements

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

        if (
            packets_mode == "virtual"
            and self.packet_data[packets_mode]["packets_df"] is None
        ):
            raise ValueError(
                "SDECPlotter doesn't have any data for virtual packets population and SDEC "
                "plot for the same was requested. Either set virtual_packet_logging: True "
                "in your configuration file to generate SDEC plot with virtual packets, or "
                "pass packets_mode='real' in your function call to generate SDEC plot with "
                "real packets."
            )

        # Store the plottable range of each spectrum property which is
        # same as entire range, initially
        self.plot_frequency_bins = self.spectrum[packets_mode][
            "spectrum_frequency_bins"
        ]
        self.plot_wavelength = self.spectrum[packets_mode][
            "spectrum_wavelength"
        ]
        self.plot_frequency = self.spectrum[packets_mode][
            "spectrum_frequency_bins"
        ][:-1]

        self.packet_wvl_range_mask = np.ones(
            self.plot_wavelength.size, dtype=bool
        )  # default value
        # Filter their plottable range based on packet_wvl_range specified
        if packet_wvl_range is not None:
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

        # Make sure number of bin edges are always one more than wavelengths
        assert self.plot_frequency_bins.size == self.plot_wavelength.size + 1

        # Calculate the area term to convert luminosity to flux
        self.lum_to_flux = 1  # default to 1 if distance is none so that this term will have no effect
        if distance is not None:
            if distance <= 0:
                raise ValueError(
                    "distance passed must be greater than 0. If you intended "
                    "to plot luminosities instead of flux, set distance=None "
                    "or don't specify distance parameter in the function call."
                )
            self.lum_to_flux = 4.0 * np.pi * (distance.to("cm")) ** 2

        # Calculate luminosities to be shown in plot
        (
            self.emission_luminosities_df,
            self.emission_species,
        ) = self._calculate_emission_luminosities(
            packets_mode=packets_mode, packet_wvl_range=packet_wvl_range
        )
        (
            self.absorption_luminosities_df,
            self.absorption_species,
        ) = self._calculate_absorption_luminosities(
            packets_mode=packets_mode, packet_wvl_range=packet_wvl_range
        )

        # Calculate the total contribution of elements
        # by summing absorption and emission
        # Only care about elements, so drop no interaction and electron scattering
        # contributions from the emitted luminosities
        self.total_luminosities_df = (
            self.absorption_luminosities_df
            + self.emission_luminosities_df.drop(["noint", "escatter"], axis=1)
        )

        # Sort the element list based on the total contribution
        sorted_list = self.total_luminosities_df.sum().sort_values(
            ascending=False
        )

        if nelements is None and self._species_list is None:
            self.species = np.array(list(self.total_luminosities_df.columns))
        elif self._species_list is not None:
            sorted_keys = list(sorted_list.keys())
            keys_to_keep = [
                key for key in sorted_keys if key in self._species_list
            ]

            df_map = {
                "total_luminosities_df": keys_to_keep,
                "emission_luminosities_df": keys_to_keep
                + ["noint", "escatter"],
                "absorption_luminosities_df": keys_to_keep,
            }

            for df_name, keys in df_map.items():
                current_df = getattr(self, df_name)
                columns_to_exclude = [
                    col for col in current_df.columns if col not in keys
                ]
                other_pos = 2 if df_name != "total_luminosities_df" else 0

                processed_df = self.process_luminosity_dataframe(
                    current_df,
                    columns_to_exclude,
                    other_column_position=other_pos,
                )
                setattr(self, df_name, processed_df)

            self.species = np.sort(self.total_luminosities_df.columns[1:])

        else:  # nelements is not None
            top_n_keys = sorted_list.keys()[:nelements]
            always_keep = ["noint", "escatter"]

            df_map = {
                "total_luminosities_df": list(top_n_keys),
                "emission_luminosities_df": list(top_n_keys) + always_keep,
                "absorption_luminosities_df": list(top_n_keys),
            }

            for df_name, keys_to_keep in df_map.items():
                current_df = getattr(self, df_name)
                columns_to_exclude = [
                    col for col in current_df.columns if col not in keys_to_keep
                ]
                other_pos = 2 if df_name != "total_luminosities_df" else 0

                processed_df = self.process_luminosity_dataframe(
                    current_df,
                    columns_to_exclude,
                    other_column_position=other_pos,
                )
                setattr(self, df_name, processed_df)

            self.species = np.sort(self.total_luminosities_df.columns[1:])

        # Final calculations
        self.photosphere_luminosity = self._calculate_photosphere_luminosity()
        self.modeled_spectrum_luminosity = (
            self.spectrum[packets_mode]["spectrum_luminosity_density_lambda"][
                self.packet_wvl_range_mask
            ]
            / self.lum_to_flux
        )

    def process_luminosity_dataframe(
        self, df, keys_to_exclude, other_column_position=0
    ):
        """
        Process a luminosity DataFrame by aggregating specified columns into an 'other' column.

        Parameters
        ----------
        df : pandas.DataFrame
            The DataFrame containing luminosity data to be processed.
        keys_to_exclude : list of str
            Column names in `df` whose data should be summed into the 'other' column and removed.
        other_column_position : int, optional
            The integer location (0-indexed) at which to insert the new 'other' column. Defaults to 0.

        Returns
        -------
        pandas.DataFrame
            A new DataFrame with the excluded columns summed into 'other' and removed from the original.
        """
        mask = np.isin(df.columns, keys_to_exclude)
        excluded_keys = df.columns[mask]

        if len(excluded_keys) > 0:
            df.insert(
                loc=other_column_position,
                column="other",
                value=df[excluded_keys].sum(axis=1),
            )
            df = df.drop(columns=excluded_keys)
        return df

    def _calculate_grouped_luminosities(
        self, packets_mode, mask, nu_column, luminosities_df
    ):
        """
        Calculate luminosities for element interactions and populate DataFrame.

        Parameters
        ----------
        packets_mode : str
            'virtual' or 'real' packets mode
        mask : np.array
            Boolean mask for wavelength filtering
        nu_column : str
            Column name containing frequency values
        luminosities_df : pd.DataFrame
            DataFrame to store results

        Returns
        -------
        tuple
            (updated luminosities_df, array of species identifiers)
        """
        # Group packets_df by atomic number of elements with which packets
        # had their last emission (interaction out)
        # or if species_list is requested then group by species id
        groupby_column = (
            "last_line_interaction_atom"
            if self._species_list is None
            else "last_line_interaction_species"
        )
        grouped = (
            self.packet_data[packets_mode]["packets_df_line_interaction"]
            .loc[mask]
            .groupby(by=groupby_column)
        )
        for identifier, group in grouped:
            weights = (
                group["energies"] / self.lum_to_flux / self.time_of_simulation
            )
            hist = np.histogram(
                group[nu_column],
                bins=self.plot_frequency_bins.value,
                weights=weights,
                density=False,
            )

            L_nu = (
                hist[0]
                * u.erg
                / u.s
                / self.spectrum[packets_mode]["spectrum_delta_frequency"]
            )
            luminosities_df[identifier] = (
                L_nu * self.plot_frequency / self.plot_wavelength
            ).value

        return luminosities_df, np.array(list(grouped.groups.keys()))

    def _calculate_luminosity_contribution(
        self, packets_mode, mask, contribution_name, luminosities_df
    ):
        """Calculate luminosity contribution for packets matching the specified mask."""
        # Histogram weights are packet luminosities or flux
        weights = (
            self.packet_data[packets_mode]["packets_df"]["energies"][
                self.packet_nu_range_mask
            ]
            / self.lum_to_flux
        ) / self.time_of_simulation
        hist = np.histogram(
            self.packet_data[packets_mode]["packets_df"]["nus"][
                self.packet_nu_range_mask
            ][mask],
            bins=self.plot_frequency_bins.value,
            weights=weights[mask],
            density=False,
        )
        L_nu = (
            hist[0]
            * u.erg
            / u.s
            / self.spectrum[packets_mode]["spectrum_delta_frequency"]
        )
        L_lambda = L_nu * self.plot_frequency / self.plot_wavelength
        luminosities_df[contribution_name] = L_lambda.value

    def _plot_traces(self, df, group_name, predefined_traces, invert_y=False):
        """Generic helper to plot traces."""
        # By specifying a common stackgroup, plotly will itself add up
        # luminosities, in order, to created stacked area chart
        base_kwargs = {
            "x": df.index,
            "stackgroup": group_name,
            **self._default_scatter_kwargs,
        }

        # Plot predefined traces first
        for colname, trace_info in predefined_traces.items():
            if colname in df.columns:
                y_data = df[colname] * (-1 if invert_y else 1)
                self.fig.add_trace(
                    go.Scatter(
                        y=y_data,
                        **trace_info,
                        **base_kwargs,
                    )
                )

        # Plot species-specific traces
        for (species_counter, identifier), species_name in zip(
            enumerate(self.species), self._species_name
        ):
            try:
                y_data = df[identifier] * (-1 if invert_y else 1)
                self.fig.add_trace(
                    go.Scatter(
                        y=y_data,
                        name=f"{species_name} {'Absorption' if invert_y else 'Emission'}",
                        fillcolor=pu.to_rgb255_string(
                            self._color_list[species_counter]
                        ),
                        hoverlabel={"namelength": -1},
                        **base_kwargs,
                    )
                )
            except KeyError:
                self._log_missing_species(
                    identifier, "absorbed" if invert_y else "emitted"
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
        self.packet_nu_range_mask = pu.create_wavelength_mask(
            self.packet_data,
            packets_mode,
            packet_wvl_range,
            df_key="packets_df",
            column_name="nus",
        )
        self.packet_nu_line_range_mask = pu.create_wavelength_mask(
            self.packet_data,
            packets_mode,
            packet_wvl_range,
            df_key="packets_df_line_interaction",
            column_name="nus",
        )

        luminosities_df = pd.DataFrame(index=self.plot_wavelength)

        # Contribution of packets which experienced no interaction
        # Mask to select packets with no interaction
        mask_noint = (
            self.packet_data[packets_mode]["packets_df"][
                "last_interaction_type"
            ][self.packet_nu_range_mask]
            == -1
        )
        self._calculate_luminosity_contribution(
            packets_mode, mask_noint, "noint", luminosities_df
        )

        # Contribution of packets which only experienced electron scattering ---
        mask_escatter = (
            self.packet_data[packets_mode]["packets_df"][
                "last_interaction_type"
            ][self.packet_nu_range_mask]
            == 1
        ) & (
            self.packet_data[packets_mode]["packets_df"][
                "last_line_interaction_in_id"
            ][self.packet_nu_range_mask]
            == -1
        )
        self._calculate_luminosity_contribution(
            packets_mode, mask_escatter, "escatter", luminosities_df
        )

        return self._calculate_grouped_luminosities(
            packets_mode=packets_mode,
            mask=self.packet_nu_line_range_mask,
            nu_column="nus",
            luminosities_df=luminosities_df,
        )

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
        self.packet_nu_line_range_mask = pu.create_wavelength_mask(
            self.packet_data,
            packets_mode,
            packet_wvl_range,
            df_key="packets_df_line_interaction",
            column_name="last_line_interaction_in_nu",
        )

        luminosities_df = pd.DataFrame(index=self.plot_wavelength)

        return self._calculate_grouped_luminosities(
            packets_mode=packets_mode,
            mask=self.packet_nu_line_range_mask,
            nu_column="last_line_interaction_in_nu",
            luminosities_df=luminosities_df,
        )

    def _calculate_photosphere_luminosity(self):
        """
        Calculate blackbody luminosity of the photosphere.

        Returns
        -------
        astropy.Quantity
            Luminosity density lambda (or Flux) of photosphere (inner boundary
            of TARDIS simulation)
        """
        bb_lam = BlackBody(
            self.t_inner,
            scale=1.0 * u.erg / (u.cm**2 * u.AA * u.s * u.sr),
        )

        L_lambda_ph = (
            bb_lam(self.plot_wavelength)
            * 4
            * np.pi**2
            * self.r_inner[0] ** 2
            * u.sr
        ).to("erg / (AA s)")

        return L_lambda_ph / self.lum_to_flux

    def generate_plot_mpl(
        self,
        packets_mode="virtual",
        packet_wvl_range=None,
        distance=None,
        observed_spectrum=None,
        show_modeled_spectrum=True,
        ax=None,
        figsize=(12, 7),
        cmapname="jet",
        nelements=None,
        species_list=None,
        blackbody_photosphere=True,
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
            Distance used to calculate flux instead of luminosity in the plot.
            It should have a length unit like m, Mpc, etc. Default value is None
        observed_spectrum : tuple or list of astropy.Quantity, optional
            Option to plot an observed spectrum in the SDEC plot. If given, the first element
            should be the wavelength and the second element should be flux,
            i.e. (wavelength, flux). The assumed units for wavelength and flux are
            angstroms and erg/(angstroms * s * cm^2), respectively. Default value is None.
        show_modeled_spectrum : bool, optional
            Whether to show modeled spectrum in SDEC Plot. Default value is
            True
        ax : matplotlib.axes._subplots.AxesSubplot or None, optional
            Axis on which to create plot. Default value is None which will
            create plot on a new figure's axis.
        figsize : tuple, optional
            Size of the matplotlib figure to display. Default value is (12, 7)
        cmapname : str, optional
            Name of matplotlib colormap to be used for showing elements.
            Default value is "jet"
        nelements: int
            Number of elements to include in plot. Determined by the
            largest contribution to total luminosity absorbed and emitted.
            Other elements are shown in silver. Default value is
            None, which displays all elements
        species_list: list of strings or None
            list of strings containing the names of species that should be included in the SDEC plots.
            Must be given in Roman numeral format. Can include specific ions, a range of ions,
            individual elements, or any combination of these:
            e.g. ['Si II', 'Ca II', 'C', 'Fe I-V']
        blackbody_photosphere: bool
            Whether to include the blackbody photosphere in the plot. Default value is True

        Returns
        -------
        matplotlib.axes._subplots.AxesSubplot
            Axis on which SDEC Plot is created
        """
        # If species_list and nelements requested, tell user that nelements is ignored
        if species_list is not None and nelements is not None:
            logger.info(
                "Both nelements and species_list were requested. Species_list takes priority; nelements is ignored"
            )

        # Parse the requested species list
        self._parse_species_list(species_list=species_list)

        # Calculate data attributes required for plotting
        # and save them in instance itself
        self._calculate_plotting_data(
            packets_mode=packets_mode,
            packet_wvl_range=packet_wvl_range,
            distance=distance,
            nelements=nelements,
        )

        if ax is None:
            self.ax = plt.figure(figsize=figsize).add_subplot(111)
        else:
            self.ax = ax

        # Get the labels in the color bar. This determines the number of unique colors
        self._make_colorbar_labels()
        # Set colormap to be used in elements of emission and absorption plots
        self.cmap = plt.get_cmap(cmapname, len(self._species_name))
        # Get the number of unqie colors
        self._make_colorbar_colors()
        self._show_colorbar_mpl()

        # Plot emission and absorption components
        self._plot_emission_mpl()
        self._plot_absorption_mpl()

        # Plot modeled spectrum
        if show_modeled_spectrum:
            self.ax.plot(
                self.plot_wavelength.value,
                self.modeled_spectrum_luminosity.value,
                "--b",
                label=f"{packets_mode.capitalize()} Spectrum",
                linewidth=1,
            )

        # Plot observed spectrum
        if observed_spectrum:
            if distance is None:
                raise ValueError(
                    """
                    Distance must be specified if an observed_spectrum is given
                    so that the model spectrum can be converted into flux space correctly.
                    """
                )

            # Convert to wavelength and luminosity units
            observed_spectrum_wavelength = observed_spectrum[0].to(u.AA)
            observed_spectrum_flux = observed_spectrum[1].to("erg/(s cm**2 AA)")

            self.ax.plot(
                observed_spectrum_wavelength.value,
                observed_spectrum_flux.value,
                "-k",
                label="Observed Spectrum",
                linewidth=1,
            )

        # Plot photosphere
        if blackbody_photosphere:
            self.ax.plot(
                self.plot_wavelength.value,
                self.photosphere_luminosity.value,
                "--r",
                label="Blackbody Photosphere",
            )

        # Set legends and labels
        xlabel = pu.axis_label_in_latex("Wavelength", u.AA)
        if distance is not None:  # Set y-axis label for flux
            ylabel = pu.axis_label_in_latex(
                "F_{\\lambda}", u.Unit("erg/(s cm**2 AA)"), only_text=False
            )
        else:  # Set y-axis label for luminosity
            ylabel = pu.axis_label_in_latex(
                "L_{\\lambda}", u.Unit("erg/(s AA)"), only_text=False
            )
        self.ax.legend(fontsize=12)
        self.ax.set_xlabel(xlabel, fontsize=12)
        self.ax.set_ylabel(
            ylabel,
            fontsize=12,
        )

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
            self.plot_wavelength.value,
            lower_level,
            upper_level,
            color="#4C4C4C",
            label="No interaction",
        )

        lower_level = upper_level
        upper_level = (
            lower_level + self.emission_luminosities_df.escatter.to_numpy()
        )

        self.ax.fill_between(
            self.plot_wavelength.value,
            lower_level,
            upper_level,
            color="#8F8F8F",
            label="Electron Scatter Only",
        )

        # If the 'other' column exists then plot it as silver
        if "other" in self.emission_luminosities_df.keys():
            lower_level = upper_level
            upper_level = (
                lower_level + self.emission_luminosities_df.other.to_numpy()
            )

            self.ax.fill_between(
                self.plot_wavelength.value,
                lower_level,
                upper_level,
                color="#C2C2C2",
                label="Other elements",
            )

        # Contribution from each element
        for species_counter, identifier in enumerate(self.species):
            try:
                lower_level = upper_level
                upper_level = (
                    lower_level
                    + self.emission_luminosities_df[identifier].to_numpy()
                )

                self.ax.fill_between(
                    self.plot_wavelength.value,
                    lower_level,
                    upper_level,
                    color=self._color_list[species_counter],
                    cmap=self.cmap,
                    linewidth=0,
                )
            except KeyError:
                # Add notifications that this species was not in the emission df
                self._log_missing_species(identifier, "emitted")

    def _plot_absorption_mpl(self):
        """Plot absorption part of the SDEC Plot using matplotlib."""
        lower_level = np.zeros(self.absorption_luminosities_df.shape[0])

        # To plot absorption part along -ve X-axis, we will start with
        # zero upper level and keep subtracting luminosities to it (lower
        # level) - fill from upper to lower level
        # If the 'other' column exists then plot it as silver
        if "other" in self.absorption_luminosities_df.keys():
            upper_level = lower_level
            lower_level = (
                upper_level - self.absorption_luminosities_df.other.to_numpy()
            )

            self.ax.fill_between(
                self.plot_wavelength.value,
                upper_level,
                lower_level,
                color="silver",
            )

        for species_counter, identifier in enumerate(self.species):
            try:
                upper_level = lower_level
                lower_level = (
                    upper_level
                    - self.absorption_luminosities_df[identifier].to_numpy()
                )

                self.ax.fill_between(
                    self.plot_wavelength.value,
                    upper_level,
                    lower_level,
                    color=self._color_list[species_counter],
                    cmap=self.cmap,
                    linewidth=0,
                )

            except KeyError:
                # Add notifications that this species was not in the emission df
                self._log_missing_species(identifier, "absorbed")

    def _show_colorbar_mpl(self):
        """Show matplotlib colorbar with labels of elements mapped to colors."""
        color_values = [
            self.cmap(species_counter / len(self._species_name))
            for species_counter in range(len(self._species_name))
        ]

        custcmap = clr.ListedColormap(color_values)
        norm = clr.Normalize(vmin=0, vmax=len(self._species_name))
        mappable = cm.ScalarMappable(norm=norm, cmap=custcmap)
        mappable.set_array(np.linspace(1, len(self._species_name) + 1, 256))
        cbar = plt.colorbar(mappable, ax=self.ax)

        bounds = np.arange(len(self._species_name)) + 0.5
        cbar.set_ticks(bounds)

        cbar.set_ticklabels(self._species_name)

    def _make_colorbar_labels(self):
        """Get the labels for the species in the colorbar."""
        if self._species_list is None:
            # If species_list is none then the labels are just elements
            species_name = [
                atomic_number2element_symbol(atomic_num)
                for atomic_num in self.species
            ]
        else:
            species_name = []
            for element in self.species:
                # Go through each species requested
                atomic_number, ion_number = divmod(element, 100)

                ion_numeral = int_to_roman(ion_number + 1)
                atomic_symbol = atomic_number2element_symbol(atomic_number)

                # if the element was requested, and not a specific ion, then
                # add the element symbol to the label list
                if (atomic_number in self._keep_colour) and (
                    atomic_symbol not in species_name
                ):
                    # compiling the label, and adding it to the list
                    label = atomic_symbol
                    species_name.append(label)
                elif atomic_number not in self._keep_colour:
                    # otherwise add the ion to the label list
                    label = f"{atomic_symbol} {ion_numeral}"
                    species_name.append(label)
        self._species_name = species_name

    def _make_colorbar_colors(self):
        """Get the colours for the species to be plotted."""
        color_list = []
        # - For elements in self._keep_colour, all ionization states share the same color
        #   (e.g., Si I, Si II, Si III all get the same color if Si's atomic number is in self._keep_colour)
        # - For elements not in self._keep_colour, each ionization state gets a new color
        for i, identifier in enumerate(self.species):
            if self._species_list is not None:
                color_counter = 0
                atomic_number = identifier // 100
                # For any element after the first one
                if i > 0:
                    previous_atomic_number = self.species[i - 1] // 100
                    # Increment color when:
                    # 1. There is a new element, OR
                    # 2. The previous element isn't in the keep_colour list
                    if (
                        previous_atomic_number != atomic_number
                        or previous_atomic_number not in self._keep_colour
                    ):
                        color_counter += 1

                color = self.cmap(color_counter / len(self._species_name))
            else:
                color = self.cmap(i / len(self.species))
            color_list.append(color)
        self._color_list = color_list

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
        nelements=None,
        species_list=None,
        blackbody_photosphere=True,
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
            Distance used to calculate flux instead of luminosity in the plot.
            It should have a length unit like m, Mpc, etc. Default value is None
        observed_spectrum : tuple or list of astropy.Quantity, optional
            Option to plot an observed spectrum in the SDEC plot. If given, the first element
            should be the wavelength and the second element should be flux,
            i.e. (wavelength, flux). The assumed units for wavelength and flux are
            angstroms and erg/(angstroms * s * cm^2), respectively. Default value is None.
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
        nelements: int
            Number of elements to include in plot. Determined by the
            largest contribution to total luminosity absorbed and emitted.
            Other elements are shown in silver. Default value is
            None, which displays all elements
        species_list: list of strings or None
            list of strings containing the names of species that should be included in the SDEC plots.
            Must be given in Roman numeral format. Can include specific ions, a range of ions,
            individual elements, or any combination of these:
            e.g. ['Si II', 'Ca II', 'C', 'Fe I-V']
        blackbody_photosphere: bool
            Whether to include the blackbody photosphere in the plot. Default value is True

        Returns
        -------
        plotly.graph_objs._figure.Figure
            Figure object on which SDEC Plot is created
        """
        # If species_list and nelements requested, tell user that nelements is ignored
        if species_list is not None and nelements is not None:
            logger.info(
                "Both nelements and species_list were requested. Species_list takes priority; nelements is ignored"
            )
        hover_props = {
            "hoverlabel": {"namelength": -1},
            "hovertemplate": "(%{x:.2f}, %{y:.3g})",
        }
        # Parse the requested species list
        self._parse_species_list(species_list=species_list)

        # Calculate data attributes required for plotting
        # and save them in instance itself
        self._calculate_plotting_data(
            packets_mode=packets_mode,
            packet_wvl_range=packet_wvl_range,
            distance=distance,
            nelements=nelements,
        )

        if fig is None:
            self.fig = go.Figure()
        else:
            self.fig = fig

        # Get the labels in the color bar. This determines the number of unique colors
        self._make_colorbar_labels()
        # Set colormap to be used in elements of emission and absorption plots
        self.cmap = plt.get_cmap(cmapname, len(self._species_name))
        # Get the number of unique colors
        self._make_colorbar_colors()

        # Plot absorption and emission components
        self._plot_emission_ply()
        self._plot_absorption_ply()

        # Plot modeled spectrum
        if show_modeled_spectrum:
            self.fig.add_trace(
                go.Scatter(
                    x=self.plot_wavelength.value,
                    y=self.modeled_spectrum_luminosity.value,
                    mode="lines",
                    line={
                        "color": "blue",
                        "width": 1,
                    },
                    name=f"{packets_mode.capitalize()} Spectrum",
                    **hover_props,
                )
            )

        # Plot observed spectrum
        if observed_spectrum:
            if distance is None:
                raise ValueError(
                    """
                    Distance must be specified if an observed_spectrum is given
                    so that the model spectrum can be converted into flux space correctly.
                    """
                )

            # Convert to wavelength and luminosity units
            observed_spectrum_wavelength = observed_spectrum[0].to(u.AA)
            observed_spectrum_flux = observed_spectrum[1].to("erg/(s cm**2 AA)")

            self.fig.add_scatter(
                x=observed_spectrum_wavelength.value,
                y=observed_spectrum_flux.value,
                name="Observed Spectrum",
                line={"color": "black", "width": 1.2},
                **hover_props,
            )

        # Plot photosphere
        if blackbody_photosphere:
            self.fig.add_trace(
                go.Scatter(
                    x=self.plot_wavelength.value,
                    y=self.photosphere_luminosity.value,
                    mode="lines",
                    line={"width": 1.5, "color": "red", "dash": "dash"},
                    name="Blackbody Photosphere",
                    **hover_props,
                )
            )

        self._show_colorbar_ply()

        # Set label and other layout options
        xlabel = pu.axis_label_in_latex("Wavelength", u.AA)
        if distance is not None:  # Set y-axis label for flux
            ylabel = pu.axis_label_in_latex(
                "F_{\\lambda}", u.Unit("erg/(s cm**2 AA)"), only_text=False
            )
        else:  # Set y-axis label for luminosity
            ylabel = pu.axis_label_in_latex(
                "L_{\\lambda}", u.Unit("erg/(s AA)"), only_text=False
            )
        self.fig.update_layout(
            xaxis={
                "title": xlabel,
                "exponentformat": "none",
            },
            yaxis={"title": ylabel, "exponentformat": "e"},
            height=graph_height,
        )

        return self.fig

    def _plot_emission_ply(self):
        """Plot emission part of the SDEC Plot using plotly."""
        self._plot_traces(
            df=self.emission_luminosities_df,
            group_name="emission",
            predefined_traces=self._predefined_traces["emission"],
            invert_y=False,
        )

    def _plot_absorption_ply(self):
        """Plot absorption part of the SDEC Plot using plotly."""
        self._plot_traces(
            df=self.absorption_luminosities_df,
            group_name="absorption",
            predefined_traces=self._predefined_traces["absorption"],
            invert_y=True,
        )

    def _show_colorbar_ply(self):
        """Show plotly colorbar with labels of elements mapped to colors."""
        # Interpolate [0, 1] range to create bins equal to number of elements
        colorscale_bins = np.linspace(0, 1, num=len(self._species_name) + 1)

        # Create a categorical colorscale [a list of (reference point, color)]
        # by mapping same reference points (excluding 1st and last bin edge)
        # twice in a row (https://plotly.com/python/colorscales/#constructing-a-discrete-or-discontinuous-color-scale)
        categorical_colorscale = []
        for species_counter in range(len(self._species_name)):
            color = pu.to_rgb255_string(
                self.cmap(colorscale_bins[species_counter])
            )
            categorical_colorscale.append(
                (colorscale_bins[species_counter], color)
            )
            categorical_colorscale.append(
                (colorscale_bins[species_counter + 1], color)
            )

        coloraxis_options = {
            "colorscale": categorical_colorscale,
            "showscale": True,
            "cmin": 0,
            "cmax": len(self._species_name),
            "colorbar": {
                "title": "Elements",
                "tickvals": np.arange(0, len(self._species_name)) + 0.5,
                "ticktext": self._species_name,
                # to change length and position of colorbar
                "len": 0.75,
                "yanchor": "top",
                "y": 0.75,
            },
        }

        # Plot an invisible one point scatter trace, to make colorbar show up
        scatter_point_idx = pu.get_mid_point_idx(self.plot_wavelength)
        self.fig.add_trace(
            go.Scatter(
                x=[self.plot_wavelength[scatter_point_idx].value],
                y=[0],
                mode="markers",
                name="Colorbar",
                showlegend=False,
                hoverinfo="skip",
                marker=dict(color=[0], opacity=0, **coloraxis_options),
            )
        )

    def _log_missing_species(self, identifier, is_absorption):
        """
        Log an informational message when a species is missing from interaction data.

        Parameters
        ----------
        identifier : int
            Species identifier, atomic number or combined atomic and ion number (Z*100+ion).
        is_absorption : bool
            True if checking absorption species, False for emission species.
        """
        interaction_type = "absorbed" if is_absorption else "emitted"
        if self._species_list is None:
            info_msg = (
                f"{atomic_number2element_symbol(identifier)}"
                f" is not in the {interaction_type} packets; skipping"
            )
        else:
            # Get the ion number and atomic number for each species
            atomic_number, ion_number = divmod(
                identifier, 100
            )  # (quotient, remainder)
            info_msg = (
                f"{atomic_number2element_symbol(atomic_number)}"
                f"{int_to_roman(ion_number + 1)}"
                f" is not in the {interaction_type} packets; skipping"
            )
        logger.info(info_msg)
