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
    element_symbol2atomic_number,
    int_to_roman,
    roman_to_int,
    species_string_to_tuple,
)
from tardis.visualization import plot_util as pu

logger = logging.getLogger(__name__)


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
        transport_state = sim.transport.transport_state
        r_inner = sim.simulation_state.geometry.r_inner_active
        t_inner = sim.simulation_state.packet_source.temperature
        time_of_simulation = (
            transport_state.packet_collection.time_of_simulation * u.s
        )

        if packets_mode == "virtual":
            return cls(
                last_interaction_type=transport_state.vpacket_tracker.last_interaction_type,
                last_line_interaction_in_id=transport_state.vpacket_tracker.last_interaction_in_id,
                last_line_interaction_out_id=transport_state.vpacket_tracker.last_interaction_out_id,
                last_line_interaction_in_nu=transport_state.vpacket_tracker.last_interaction_in_nu,
                lines_df=lines_df,
                packet_nus=u.Quantity(
                    transport_state.vpacket_tracker.nus, "Hz"
                ),
                packet_energies=u.Quantity(
                    transport_state.vpacket_tracker.energies, "erg"
                ),
                r_inner=r_inner,
                spectrum_delta_frequency=transport_state.spectrum_virtual.delta_frequency,
                spectrum_frequency_bins=transport_state.spectrum_virtual._frequency,
                spectrum_luminosity_density_lambda=transport_state.spectrum_virtual.luminosity_density_lambda,
                spectrum_wavelength=transport_state.spectrum_virtual.wavelength,
                t_inner=t_inner,
                time_of_simulation=time_of_simulation,
            )

        elif packets_mode == "real":
            # Packets-specific properties need to be only for those packets
            # which got emitted
            transport_state = sim.transport.transport_state
            return cls(
                last_interaction_type=transport_state.last_interaction_type[
                    transport_state.emitted_packet_mask
                ],
                last_line_interaction_in_id=transport_state.last_line_interaction_in_id[
                    transport_state.emitted_packet_mask
                ],
                last_line_interaction_out_id=transport_state.last_line_interaction_out_id[
                    transport_state.emitted_packet_mask
                ],
                last_line_interaction_in_nu=transport_state.last_interaction_in_nu[
                    transport_state.emitted_packet_mask
                ],
                lines_df=lines_df,
                packet_nus=transport_state.packet_collection.output_nus[
                    transport_state.emitted_packet_mask
                ],
                packet_energies=transport_state.packet_collection.output_energies[
                    transport_state.emitted_packet_mask
                ],
                r_inner=r_inner,
                spectrum_delta_frequency=transport_state.spectrum.delta_frequency,
                spectrum_frequency_bins=transport_state.spectrum._frequency,
                spectrum_luminosity_density_lambda=transport_state.spectrum.luminosity_density_lambda,
                spectrum_wavelength=transport_state.spectrum.wavelength,
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
                hdf["/simulation/simulation_state/r_inner"].to_numpy(), "cm"
            )  # Convert pd.Series to np.array to construct quantity from it
            t_inner = u.Quantity(
                hdf["/simulation/simulation_state/scalars"].t_inner, "K"
            )
            time_of_simulation = u.Quantity(
                hdf["/simulation/transport/scalars"].time_of_simulation, "s"
            )

            if packets_mode == "virtual":
                return cls(
                    last_interaction_type=hdf[
                        "/simulation/transport/transport_state/virt_packet_last_interaction_type"
                    ],
                    last_line_interaction_in_id=hdf[
                        "/simulation/transport/transport_state/virt_packet_last_line_interaction_in_id"
                    ],
                    last_line_interaction_out_id=hdf[
                        "/simulation/transport/transport_state/virt_packet_last_line_interaction_out_id"
                    ],
                    last_line_interaction_in_nu=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/virt_packet_last_interaction_in_nu"
                        ].to_numpy(),
                        "Hz",
                    ),
                    lines_df=lines_df,
                    packet_nus=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/virt_packet_nus"
                        ].to_numpy(),
                        "Hz",
                    ),
                    packet_energies=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/virt_packet_energies"
                        ].to_numpy(),
                        "erg",
                    ),
                    r_inner=r_inner,
                    spectrum_delta_frequency=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/spectrum_virtual/scalars"
                        ].delta_frequency,
                        "Hz",
                    ),
                    spectrum_frequency_bins=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/spectrum_virtual/_frequency"
                        ].to_numpy(),
                        "Hz",
                    ),
                    spectrum_luminosity_density_lambda=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/spectrum_virtual/luminosity_density_lambda"
                        ].to_numpy(),
                        "erg / s cm",  # luminosity_density_lambda is saved in hdf in CGS
                    ).to("erg / s AA"),
                    spectrum_wavelength=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/spectrum_virtual/wavelength"
                        ].to_numpy(),
                        "cm",  # wavelength is saved in hdf in CGS
                    ).to("AA"),
                    t_inner=t_inner,
                    time_of_simulation=time_of_simulation,
                )

            elif packets_mode == "real":
                emitted_packet_mask = hdf[
                    "/simulation/transport/transport_state/emitted_packet_mask"
                ].to_numpy()
                return cls(
                    # First convert series read from hdf to array before masking
                    # to eliminate index info which creates problems otherwise
                    last_interaction_type=hdf[
                        "/simulation/transport/transport_state/last_interaction_type"
                    ].to_numpy()[emitted_packet_mask],
                    last_line_interaction_in_id=hdf[
                        "/simulation/transport/transport_state/last_line_interaction_in_id"
                    ].to_numpy()[emitted_packet_mask],
                    last_line_interaction_out_id=hdf[
                        "/simulation/transport/transport_state/last_line_interaction_out_id"
                    ].to_numpy()[emitted_packet_mask],
                    last_line_interaction_in_nu=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/last_interaction_in_nu"
                        ].to_numpy()[emitted_packet_mask],
                        "Hz",
                    ),
                    lines_df=lines_df,
                    packet_nus=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/output_nu"
                        ].to_numpy()[emitted_packet_mask],
                        "Hz",
                    ),
                    packet_energies=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/output_energy"
                        ].to_numpy()[emitted_packet_mask],
                        "erg",
                    ),
                    r_inner=r_inner,
                    spectrum_delta_frequency=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/spectrum/scalars"
                        ].delta_frequency,
                        "Hz",
                    ),
                    spectrum_frequency_bins=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/spectrum/_frequency"
                        ].to_numpy(),
                        "Hz",
                    ),
                    spectrum_luminosity_density_lambda=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/spectrum/luminosity_density_lambda"
                        ].to_numpy(),
                        "erg / s cm",
                    ).to("erg / s AA"),
                    spectrum_wavelength=u.Quantity(
                        hdf[
                            "/simulation/transport/transport_state/spectrum/wavelength"
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
        if sim.transport.transport_state.virt_logging:
            return cls(
                {
                    "virtual": SDECData.from_simulation(sim, "virtual"),
                    "real": SDECData.from_simulation(sim, "real"),
                }
            )
        else:
            return cls(
                {
                    "virtual": None,
                    "real": SDECData.from_simulation(sim, "real"),
                }
            )

    @classmethod
    def from_hdf(cls, hdf_fpath, packets_mode=None):
        """
        Create an instance of SDECPlotter from a simulation HDF file.

        Parameters
        ----------
        hdf_fpath : str
            Valid path to the HDF file where simulation is saved
        packets_mode : {'virtual', 'real'}, optional
            Mode of packets to be considered, either real or virtual. If not
            specified, both modes are returned

        Returns
        -------
        SDECPlotter
        """
        assert packets_mode in [None, "virtual", "real"], (
            "Invalid value passed to packets_mode. Only "
            "allowed values are 'virtual', 'real' or None"
        )
        if packets_mode == "virtual":
            return cls(
                {
                    "virtual": SDECData.from_hdf(hdf_fpath, "virtual"),
                    "real": None,
                }
            )
        elif packets_mode == "real":
            return cls(
                {
                    "virtual": None,
                    "real": SDECData.from_hdf(hdf_fpath, "real"),
                }
            )
        else:
            return cls(
                {
                    "virtual": SDECData.from_hdf(hdf_fpath, "virtual"),
                    "real": SDECData.from_hdf(hdf_fpath, "real"),
                }
            )

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
            # check if there are any digits in the species list. If there are, then exit.
            # species_list should only contain species in the Roman numeral
            # format, e.g. Si II, and each ion must contain a space
            if any(char.isdigit() for char in " ".join(species_list)) is True:
                raise ValueError(
                    "All species must be in Roman numeral form, e.g. Si II"
                )
            else:
                full_species_list = []
                for species in species_list:
                    # check if a hyphen is present. If it is, then it indicates a
                    # range of ions. Add each ion in that range to the list as a new entry
                    if "-" in species:
                        # split the string on spaces. First thing in the list is then the element
                        element = species.split(" ")[0]
                        # Next thing is the ion range
                        # convert the requested ions into numerals
                        first_ion_numeral = roman_to_int(
                            species.split(" ")[-1].split("-")[0]
                        )
                        second_ion_numeral = roman_to_int(
                            species.split(" ")[-1].split("-")[-1]
                        )
                        # add each ion between the two requested into the species list
                        for ion_number in np.arange(
                            first_ion_numeral, second_ion_numeral + 1
                        ):
                            full_species_list.append(
                                f"{element} {int_to_roman(ion_number)}"
                            )
                    else:
                        # Otherwise it's either an element or ion so just add to the list
                        full_species_list.append(species)

                # full_species_list is now a list containing each individual species requested
                # e.g. it parses species_list = [Si I - V] into species_list = [Si I, Si II, Si III, Si IV, Si V]
                self._full_species_list = full_species_list
                requested_species_ids = []
                keep_colour = []

                # go through each of the requested species. Check whether it is
                # an element or ion (ions have spaces). If it is an element,
                # add all possible ions to the ions list. Otherwise just add
                # the requested ion
                for species in full_species_list:
                    if " " in species:
                        requested_species_ids.append(
                            [
                                species_string_to_tuple(species)[0] * 100
                                + species_string_to_tuple(species)[1]
                            ]
                        )
                    else:
                        atomic_number = element_symbol2atomic_number(species)
                        requested_species_ids.append(
                            [
                                atomic_number * 100 + ion_number
                                for ion_number in np.arange(atomic_number)
                            ]
                        )
                        # add the atomic number to a list so you know that this element should
                        # have all species in the same colour, i.e. it was requested like
                        # species_list = [Si]
                        keep_colour.append(atomic_number)
                requested_species_ids = [
                    species_id
                    for temp_list in requested_species_ids
                    for species_id in temp_list
                ]

                self._species_list = requested_species_ids
                self._keep_colour = keep_colour
        else:
            self._species_list = None

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

        if packets_mode == "virtual" and self.data[packets_mode] is None:
            raise ValueError(
                "SDECPlotter doesn't have any data for virtual packets population and SDEC "
                "plot for the same was requested. Either set virtual_packet_logging: True "
                "in your configuration file to generate SDEC plot with virtual packets, or "
                "pass packets_mode='real' in your function call to generate SDEC plot with "
                "real packets."
            )

        # Store the plottable range of each spectrum property which is
        # same as entire range, initially
        self.plot_frequency_bins = self.data[
            packets_mode
        ].spectrum_frequency_bins
        self.plot_wavelength = self.data[packets_mode].spectrum_wavelength
        self.plot_frequency = self.data[packets_mode].spectrum_frequency

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

        # If nelements and species_list are not included, the list of elements is just all elements
        if nelements is None and self._species_list is None:
            self.species = np.array(list(self.total_luminosities_df.keys()))
        elif self._species_list is not None:
            # Compare the species present in the model to those in the requested list
            # Mask out those that aren't in the model
            mask = np.in1d(
                np.array(list(sorted_list.keys())), self._species_list
            )
            # If species_list is included then create a new column which is the sum
            # of all other species i.e. those that aren't in the requested list
            self.total_luminosities_df.insert(
                loc=0,
                column="other",
                value=self.total_luminosities_df[sorted_list.keys()[~mask]].sum(
                    axis=1
                ),
            )
            # Then drop all of the individual columns for species included in 'other'
            self.total_luminosities_df.drop(
                sorted_list.keys()[~mask], inplace=True, axis=1
            )
            # Repeat this for the emission and absorption dfs
            # This will require creating a temporary list that includes 'noint' and 'escatter'
            # packets, because you don't want them dropped or included in 'other'
            temp = list(self._species_list)
            temp.append("noint")
            temp.append("escatter")
            mask = np.in1d(
                np.array(list(self.emission_luminosities_df.keys())), temp
            )
            # If species_list is included then create a new column which is the sum
            # of all other species i.e. those that aren't in the requested list
            self.emission_luminosities_df.insert(
                loc=0,
                column="other",
                value=self.emission_luminosities_df[
                    self.emission_luminosities_df.keys()[~mask]
                ].sum(axis=1),
            )
            # Need to add a new value to the mask array for the 'other' column just added
            mask = np.insert(mask, 0, True)
            # Then drop all of the individual columns for species included in 'other'
            self.emission_luminosities_df.drop(
                self.emission_luminosities_df.keys()[~mask],
                inplace=True,
                axis=1,
            )

            temp = list(self._species_list)
            mask = np.in1d(
                np.array(list(self.absorption_luminosities_df.keys())), temp
            )
            # If species_list is included then create a new column which is the sum
            # of all other species i.e. those that aren't in the requested list
            self.absorption_luminosities_df.insert(
                loc=0,
                column="other",
                value=self.absorption_luminosities_df[
                    self.absorption_luminosities_df.keys()[~mask]
                ].sum(axis=1),
            )
            # Need to add a new value to the mask array for the 'other' column just added
            mask = np.insert(mask, 0, True)
            # Then drop all of the individual columns for species included in 'other'
            self.absorption_luminosities_df.drop(
                self.absorption_luminosities_df.keys()[~mask],
                inplace=True,
                axis=1,
            )

            # Get the list of species in the model
            # Index from 1: to avoid the 'other' column
            self.species = np.sort(self.total_luminosities_df.keys()[1:])

        else:
            # If nelements is included then create a new column which is the sum
            # of all other elements, i.e. those that aren't in the top contributing nelements
            self.total_luminosities_df.insert(
                loc=0,
                column="other",
                value=self.total_luminosities_df[
                    sorted_list.keys()[nelements:]
                ].sum(axis=1),
            )
            # Then drop all of the individual columns for elements included in 'other'
            self.total_luminosities_df.drop(
                sorted_list.keys()[nelements:], inplace=True, axis=1
            )
            # If nelements is included then create a new column which is the sum
            # of all other elements, i.e. those that aren't in the top contributing nelements
            self.emission_luminosities_df.insert(
                loc=2,
                column="other",
                value=self.emission_luminosities_df[
                    sorted_list.keys()[nelements:]
                ].sum(axis=1),
            )
            # Then drop all of the individual columns for elements included in 'other'
            self.emission_luminosities_df.drop(
                sorted_list.keys()[nelements:], inplace=True, axis=1
            )
            # If nelements is included then create a new column which is the sum
            # of all other elements, i.e. those that aren't in the top contributing nelements
            self.absorption_luminosities_df.insert(
                loc=2,
                column="other",
                value=self.absorption_luminosities_df[
                    sorted_list.keys()[nelements:]
                ].sum(axis=1),
            )
            # Then drop all of the individual columns for elements included in 'other'
            self.absorption_luminosities_df.drop(
                sorted_list.keys()[nelements:], inplace=True, axis=1
            )
            # Get the list of species in the model
            # Index from 1: to avoid the 'other' column
            self.species = np.sort(self.total_luminosities_df.keys()[1:])

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
            bins=self.plot_frequency_bins.value,
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
            bins=self.plot_frequency_bins.value,
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
        # or if species_list is requested then group by species id
        if self._species_list is None:
            packets_df_grouped = (
                self.data[packets_mode]
                .packets_df_line_interaction.loc[self.packet_nu_line_range_mask]
                .groupby(by="last_line_interaction_atom")
            )
        else:
            packets_df_grouped = (
                self.data[packets_mode]
                .packets_df_line_interaction.loc[self.packet_nu_line_range_mask]
                .groupby(by="last_line_interaction_species")
            )

        # Contribution of each species with which packets interacted ----------
        for identifier, group in packets_df_grouped:
            # Histogram of specific species
            hist_el = np.histogram(
                group["nus"],
                bins=self.plot_frequency_bins.value,
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

            luminosities_df[identifier] = L_lambda_el.value

        # Create an array of the species with which packets interacted
        emission_species_present = np.array(
            list(packets_df_grouped.groups.keys())
        )

        return luminosities_df, emission_species_present

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
        # or if species_list is requested then group by species id
        if self._species_list is None:
            packets_df_grouped = (
                self.data[packets_mode]
                .packets_df_line_interaction.loc[self.packet_nu_line_range_mask]
                .groupby(by="last_line_interaction_atom")
            )
        else:
            packets_df_grouped = (
                self.data[packets_mode]
                .packets_df_line_interaction.loc[self.packet_nu_line_range_mask]
                .groupby(by="last_line_interaction_species")
            )

        for identifier, group in packets_df_grouped:
            # Histogram of specific species
            hist_el = np.histogram(
                group["last_line_interaction_in_nu"],
                bins=self.plot_frequency_bins.value,
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

            luminosities_df[identifier] = L_lambda_el.value

        absorption_species_present = np.array(
            list(packets_df_grouped.groups.keys())
        )

        return luminosities_df, absorption_species_present

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
        bb_lam = BlackBody(
            self.data[packets_mode].t_inner,
            scale=1.0 * u.erg / (u.cm**2 * u.AA * u.s * u.sr),
        )

        L_lambda_ph = (
            bb_lam(self.plot_wavelength)
            * 4
            * np.pi**2
            * self.data[packets_mode].r_inner[0] ** 2
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
        self.cmap = cm.get_cmap(cmapname, len(self._species_name))
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

            observed_spectrum_wavelength = None
            observed_spectrum_flux = None

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
        self.ax.legend(fontsize=12)
        self.ax.set_xlabel(r"Wavelength $[\mathrm{\AA}]$", fontsize=12)
        if distance is not None:  # Set y-axis label for flux
            self.ax.set_ylabel(
                r"$F_{\lambda}$ [erg $\mathrm{s^{-1}}$ $\mathrm{cm^{-2}}$ $\mathrm{\AA^{-1}}$]",
                fontsize=12,
            )
        else:  # Set y-axis label for luminosity
            self.ax.set_ylabel(
                r"$L_{\lambda}$ [erg $\mathrm{s^{-1}}$ $\mathrm{\AA^{-1}}$]",
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
                if self._species_list is None:
                    info_msg = (
                        f"{atomic_number2element_symbol(identifier)}"
                        f" is not in the emitted packets; skipping"
                    )
                    logger.info(info_msg)
                else:
                    # Get the ion number and atomic number for each species
                    ion_number = identifier % 100
                    atomic_number = (identifier - ion_number) / 100

                    info_msg = (
                        f"{atomic_number2element_symbol(atomic_number)}"
                        f"{int_to_roman(ion_number + 1)}"
                        f" is not in the emitted packets; skipping"
                    )
                    logger.info(info_msg)

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
                if self._species_list is None:
                    info_msg = (
                        f"{atomic_number2element_symbol(identifier)}"
                        f" is not in the absorbed packets; skipping"
                    )
                    logger.info(info_msg)
                else:
                    # Get the ion number and atomic number for each species
                    ion_number = identifier % 100
                    atomic_number = (identifier - ion_number) / 100

                    info_msg = (
                        f"{atomic_number2element_symbol(atomic_number)}"
                        f"{int_to_roman(ion_number + 1)}"
                        f" is not in the absorbed packets; skipping"
                    )
                    logger.info(info_msg)

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
            for species in self.species:
                # Go through each species requested
                ion_number = species % 100
                atomic_number = (species - ion_number) / 100

                ion_numeral = int_to_roman(ion_number + 1)
                atomic_symbol = atomic_number2element_symbol(atomic_number)

                # if the element was requested, and not a specific ion, then
                # add the element symbol to the label list
                if (atomic_number in self._keep_colour) & (
                    atomic_symbol not in species_name
                ):
                    # compiling the label, and adding it to the list
                    label = f"{atomic_symbol}"
                    species_name.append(label)
                elif atomic_number not in self._keep_colour:
                    # otherwise add the ion to the label list
                    label = f"{atomic_symbol} {ion_numeral}"
                    species_name.append(label)

        self._species_name = species_name

    def _make_colorbar_colors(self):
        """Get the colours for the species to be plotted."""
        # the colours depends on the species present in the model and what's requested
        # some species need to be shown in the same colour, so the exact colours have to be
        # worked out

        color_list = []

        # Colors for each element
        # Create new variables to keep track of the last atomic number that was plotted
        # This is used when plotting species in case an element was given in the list
        # This is to ensure that all ions of that element are grouped together
        # ii is to track the colour index
        # e.g. if Si is given in species_list, this is to ensure Si I, Si II, etc. all have the same colour
        color_counter = 0
        previous_atomic_number = 0
        for species_counter, identifier in enumerate(self.species):
            if self._species_list is not None:
                # Get the ion number and atomic number for each species
                ion_number = identifier % 100
                atomic_number = (identifier - ion_number) / 100
                if previous_atomic_number == 0:
                    # If this is the first species being plotted, then take note of the atomic number
                    # don't update the colour index
                    previous_atomic_number = atomic_number
                elif previous_atomic_number in self._keep_colour:
                    # If the atomic number is in the list of elements that should all be plotted in the same colour
                    # then don't update the colour index if this element has been plotted already
                    if previous_atomic_number == atomic_number:
                        previous_atomic_number = atomic_number
                    else:
                        # Otherwise, increase the colour counter by one, because this is a new element
                        color_counter = color_counter + 1
                        previous_atomic_number = atomic_number
                else:
                    # If this is just a normal species that was requested then increment the colour index
                    color_counter = color_counter + 1
                    previous_atomic_number = atomic_number
                # Calculate the colour of this species
                color = self.cmap(color_counter / len(self._species_name))

            else:
                # If you're not using species list then this is just a fraction based on the total
                # number of columns in the dataframe
                color = self.cmap(species_counter / len(self.species))

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
        self.cmap = cm.get_cmap(cmapname, len(self._species_name))
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
                    hovertemplate="(%{x:.2f}, %{y:.3g})",
                    hoverlabel={"namelength": -1},
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

            observed_spectrum_wavelength = None
            observed_spectrum_flux = None

            # Convert to wavelength and luminosity units
            observed_spectrum_wavelength = observed_spectrum[0].to(u.AA)
            observed_spectrum_flux = observed_spectrum[1].to("erg/(s cm**2 AA)")

            self.fig.add_scatter(
                x=observed_spectrum_wavelength.value,
                y=observed_spectrum_flux.value,
                name="Observed Spectrum",
                line={"color": "black", "width": 1.2},
                hoverlabel={"namelength": -1},
                hovertemplate="(%{x:.2f}, %{y:.3g})",
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
                    hoverlabel={"namelength": -1},
                    hovertemplate="(%{x:.2f}, %{y:.3g})",
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
                fillcolor="#4C4C4C",
                stackgroup="emission",
                hovertemplate="(%{x:.2f}, %{y:.3g})",
            )
        )

        self.fig.add_trace(
            go.Scatter(
                x=self.emission_luminosities_df.index,
                y=self.emission_luminosities_df.escatter,
                mode="none",
                name="Electron Scatter Only",
                fillcolor="#8F8F8F",
                stackgroup="emission",
                hoverlabel={"namelength": -1},
                hovertemplate="(%{x:.2f}, %{y:.3g})",
            )
        )

        # If 'other' column exists then plot as silver
        if "other" in self.emission_luminosities_df.keys():
            self.fig.add_trace(
                go.Scatter(
                    x=self.emission_luminosities_df.index,
                    y=self.emission_luminosities_df.other,
                    mode="none",
                    name="Other elements",
                    fillcolor="#C2C2C2",
                    stackgroup="emission",
                    hovertemplate="(%{x:.2f}, %{y:.3g})",
                )
            )

        # Contribution from each element
        for (species_counter, identifier), species_name in zip(
            enumerate(self.species), self._species_name
        ):
            try:
                self.fig.add_trace(
                    go.Scatter(
                        x=self.emission_luminosities_df.index,
                        y=self.emission_luminosities_df[identifier],
                        mode="none",
                        name=species_name + " Emission",
                        hovertemplate=f"<b>{species_name:s} Emission<br>"  # noqa: ISC003
                        + "(%{x:.2f}, %{y:.3g})<extra></extra>",
                        fillcolor=self.to_rgb255_string(
                            self._color_list[species_counter]
                        ),
                        stackgroup="emission",
                        showlegend=False,
                        hoverlabel={"namelength": -1},
                    )
                )
            except KeyError:
                # Add notifications that this species was not in the emission df
                if self._species_list is None:
                    info_msg = (
                        f"{atomic_number2element_symbol(identifier)}"
                        f" is not in the emitted packets; skipping"
                    )
                    logger.info(info_msg)
                else:
                    # Get the ion number and atomic number for each species
                    ion_number = identifier % 100
                    atomic_number = (identifier - ion_number) / 100

                    info_msg = (
                        f"{atomic_number2element_symbol(atomic_number)}"
                        f"{int_to_roman(ion_number + 1)}"
                        f" is not in the emitted packets; skipping"
                    )
                    logger.info(info_msg)

    def _plot_absorption_ply(self):
        """Plot absorption part of the SDEC Plot using plotly."""
        # If 'other' column exists then plot as silver
        if "other" in self.absorption_luminosities_df.keys():
            self.fig.add_trace(
                go.Scatter(
                    x=self.absorption_luminosities_df.index,
                    # to plot absorption luminosities along negative y-axis
                    y=self.absorption_luminosities_df.other * -1,
                    mode="none",
                    name="Other elements",
                    fillcolor="#C2C2C2",
                    stackgroup="absorption",
                    showlegend=False,
                    hovertemplate="(%{x:.2f}, %{y:.3g})",
                )
            )

        for (species_counter, identifier), species_name in zip(
            enumerate(self.species), self._species_name
        ):
            try:
                self.fig.add_trace(
                    go.Scatter(
                        x=self.absorption_luminosities_df.index,
                        # to plot absorption luminosities along negative y-axis
                        y=self.absorption_luminosities_df[identifier] * -1,
                        mode="none",
                        name=species_name + " Absorption",
                        hovertemplate=f"<b>{species_name:s} Absorption<br>"  # noqa: ISC003
                        + "(%{x:.2f}, %{y:.3g})<extra></extra>",
                        fillcolor=self.to_rgb255_string(
                            self._color_list[species_counter]
                        ),
                        stackgroup="absorption",
                        showlegend=False,
                        hoverlabel={"namelength": -1},
                    )
                )

            except KeyError:
                # Add notifications that this species was not in the emission df
                if self._species_list is None:
                    info_msg = (
                        f"{atomic_number2element_symbol(identifier)}"
                        f" is not in the absorbed packets; skipping"
                    )
                    logger.info(info_msg)
                else:
                    # Get the ion number and atomic number for each species
                    ion_number = identifier % 100
                    atomic_number = (identifier - ion_number) / 100

                    info_msg = (
                        f"{atomic_number2element_symbol(atomic_number)}"
                        f"{int_to_roman(ion_number + 1)}"
                        f" is not in the absorbed packets; skipping"
                    )
                    logger.info(info_msg)

    def _show_colorbar_ply(self):
        """Show plotly colorbar with labels of elements mapped to colors."""
        # Interpolate [0, 1] range to create bins equal to number of elements
        colorscale_bins = np.linspace(0, 1, num=len(self._species_name) + 1)

        # Create a categorical colorscale [a list of (reference point, color)]
        # by mapping same reference points (excluding 1st and last bin edge)
        # twice in a row (https://plotly.com/python/colorscales/#constructing-a-discrete-or-discontinuous-color-scale)
        categorical_colorscale = []
        for species_counter in range(len(self._species_name)):
            color = self.to_rgb255_string(
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
