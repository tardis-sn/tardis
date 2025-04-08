"""Utility functions to be used in plotting."""

import re

import astropy.units as u
import numpy as np
import pandas as pd

from tardis.util.base import (
    atomic_number2element_symbol,
    element_symbol2atomic_number,
    int_to_roman,
    roman_to_int,
    species_string_to_tuple,
)

from tardis.util.base import (
    atomic_number2element_symbol,
    element_symbol2atomic_number,
    int_to_roman,
    roman_to_int,
    species_string_to_tuple,
)


def axis_label_in_latex(label_text, unit, only_text=True):
    """
    Get axis label for plotly plots that can show units in latex.

    Parameters
    ----------
    label_text : str
        Text to show on label, may be expressed in latex
    unit : astropy.units
        Unit of the label which needs to be expressed in latex
    only_text : bool
        If label_text is expressed purely in text (i.e. without
        using latex) or not. Default value is True

    Returns
    -------
    str
        Latex string for label renderable by plotly
    """
    unit_in_latex = unit.to_string("latex_inline").strip("$")

    # If present, place s^{-1} just after erg
    if "erg" in unit_in_latex and "s^{-1}" in unit_in_latex:
        constituent_units = (
            re.compile(r"\\mathrm\{(.*)\}")
            .findall(unit_in_latex)[0]
            .split("\\,")
        )
        constituent_units.remove("s^{-1}")
        constituent_units.insert(constituent_units.index("erg") + 1, "s^{-1}")
        constituent_units_string = "\\,".join(constituent_units)
        unit_in_latex = f"\\mathrm{{{constituent_units_string}}}"

    if only_text:
        return f"$\\text{{{label_text}}}\\,[{unit_in_latex}]$"
    else:
        return f"${label_text}\\,[{unit_in_latex}]$"


def get_mid_point_idx(arr):
    """
    Get index of the middle point of a sorted array (ascending or descending).

    The values in array may not be evenly distributed so it picks the middle
    point not by index but by their values.

    Parameters
    ----------
    arr : np.array

    Returns
    -------
    int
    """
    mid_value = (arr[0] + arr[-1]) / 2
    return np.abs(arr - mid_value).argmin()


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


def extract_and_process_packet_data(simulation, packets_mode):
    """
    Extract and process packet data from the simulation object.

    This includes converting the packet data into a DataFrame and appending
    processed line interaction information.

    Parameters
    ----------
    simulation : tardis.simulation.BaseSimulation
        Full TARDIS simulation object containing transport state, plasma,
        and atomic data.

    packets_mode : str
        Type of packets to extract:
        - 'virtual': Use virtual packet tracker.
        - 'real': Use emitted real packets.

    Returns
    -------
    dict
        Dictionary containing raw packet data, the full DataFrame `packets_df`,
        and a filtered `packets_df_line_interaction` with line interaction info.
    """
    transport_state = simulation.transport.transport_state
    lines_df = simulation.plasma.atomic_data.lines.reset_index().set_index(
        "line_id"
    )

    if packets_mode == "virtual":
        vpacket_tracker = transport_state.vpacket_tracker
        packet_data = {
            "last_interaction_type": vpacket_tracker.last_interaction_type,
            "last_line_interaction_in_id": vpacket_tracker.last_interaction_in_id,
            "last_line_interaction_out_id": vpacket_tracker.last_interaction_out_id,
            "last_line_interaction_in_nu": vpacket_tracker.last_interaction_in_nu,
            "last_interaction_in_r": vpacket_tracker.last_interaction_in_r,
            "nus": u.Quantity(vpacket_tracker.nus, "Hz"),
            "energies": u.Quantity(vpacket_tracker.energies, "erg"),
            "lambdas": u.Quantity(vpacket_tracker.nus, "Hz").to(
                "angstrom", u.spectral()
            ),
        }
    else:
        mask = transport_state.emitted_packet_mask
        packet_nus = u.Quantity(
            transport_state.packet_collection.output_nus[mask], u.Hz
        )
        packet_data = {
            "last_interaction_type": transport_state.last_interaction_type[
                mask
            ],
            "last_line_interaction_in_id": transport_state.last_line_interaction_in_id[
                mask
            ],
            "last_line_interaction_out_id": transport_state.last_line_interaction_out_id[
                mask
            ],
            "last_line_interaction_in_nu": transport_state.last_interaction_in_nu[
                mask
            ],
            "last_interaction_in_r": transport_state.last_interaction_in_r[
                mask
            ],
            "nus": packet_nus,
            "energies": transport_state.packet_collection.output_energies[mask],
            "lambdas": packet_nus.to("angstrom", u.spectral()),
        }

    packet_data["packets_df"] = pd.DataFrame(packet_data)
    process_line_interactions(packet_data, lines_df)
    return packet_data


def process_line_interactions(packet_data, lines_df):
    """
    Add line interaction metadata to the packet DataFrame.

    Filters packets that experienced a line interaction and computes:
    - Atomic number of the last interaction (out).
    - Species ID (Z * 100 + ion number).

    Parameters
    ----------
    packet_data : dict
        Dictionary containing a 'packets_df' DataFrame.

    lines_df : pandas.DataFrame
        DataFrame with 'atomic_number' and 'ion_number' indexed by line ID.
    """
    packets_df = packet_data["packets_df"]

    if packets_df is not None:
        # Create dataframe of packets that experience line interaction
        line_mask = (packets_df["last_interaction_type"] > -1) & (
            packets_df["last_line_interaction_in_id"] > -1
        )
        packet_data["packets_df_line_interaction"] = packets_df.loc[
            line_mask
        ].copy()

        # Add columns for atomic number of last interaction out
        packet_data["packets_df_line_interaction"][
            "last_line_interaction_atom"
        ] = (
            lines_df["atomic_number"]
            .iloc[
                packet_data["packets_df_line_interaction"][
                    "last_line_interaction_out_id"
                ]
            ]
            .to_numpy()
        )

        # Add columns for the species ID of last interaction
        packet_data["packets_df_line_interaction"][
            "last_line_interaction_species"
        ] = (
            lines_df["atomic_number"]
            .iloc[
                packet_data["packets_df_line_interaction"][
                    "last_line_interaction_out_id"
                ]
            ]
            .to_numpy()
            * 100
            + lines_df["ion_number"]
            .iloc[
                packet_data["packets_df_line_interaction"][
                    "last_line_interaction_out_id"
                ]
            ]
            .to_numpy()
        )


def extract_and_process_packet_data_hdf(hdf, packets_mode):
    """
    Extract and process packet data from an HDF file.

    This function retrieves packet data from an HDF file and processes it
    based on the specified packet mode (either "virtual" or "real"). The
    extracted data is organized into a dictionary and further processed
    to include line interaction information.

    Parameters
    ----------
    hdf : h5py.File or dict-like
        The HDF file object containing the simulation data.
    packets_mode : str
        The mode of packets to process. Can be "virtual" for virtual packets
        or any other value for real packets.

    Returns
    -------
    dict
        A dictionary containing the processed packet data.

    Raises
    ------
    KeyError
        If required keys are missing in the HDF file.
    ValueError
        If an invalid `packets_mode` is provided.
    """
    lines_df = (
        hdf["/simulation/plasma/lines"].reset_index().set_index("line_id")
    )
    if packets_mode == "virtual":
        packet_prefix = "/simulation/transport/transport_state/virt_packet"
        packet_data = {
            "last_interaction_type": hdf[
                f"{packet_prefix}_last_interaction_type"
            ],
            "last_line_interaction_in_id": hdf[
                f"{packet_prefix}_last_interaction_in_id"
            ],
            "last_line_interaction_out_id": hdf[
                f"{packet_prefix}_last_interaction_out_id"
            ],
            "last_line_interaction_in_nu": u.Quantity(
                hdf[f"{packet_prefix}_last_interaction_in_nu"].to_numpy(), "Hz"
            ),
            "last_interaction_in_r": u.Quantity(
                hdf[f"{packet_prefix}_last_interaction_in_r"].to_numpy(), "cm"
            ),
            "packet_nus": u.Quantity(
                hdf[f"{packet_prefix}_nus"].to_numpy(), "Hz"
            ),
            "packet_energies": u.Quantity(
                hdf[f"{packet_prefix}_energies"].to_numpy(), "erg"
            ),
        }
    else:  # real packets
        emitted_packet_mask = hdf[
            "/simulation/transport/transport_state/emitted_packet_mask"
        ].to_numpy()
        packet_prefix = "/simulation/transport/transport_state"
        packet_data = {
            "last_interaction_type": hdf[
                f"{packet_prefix}/last_interaction_type"
            ].to_numpy()[emitted_packet_mask],
            "last_line_interaction_in_id": hdf[
                f"{packet_prefix}/last_line_interaction_in_id"
            ].to_numpy()[emitted_packet_mask],
            "last_line_interaction_out_id": hdf[
                f"{packet_prefix}/last_line_interaction_out_id"
            ].to_numpy()[emitted_packet_mask],
            "last_line_interaction_in_nu": u.Quantity(
                hdf[f"{packet_prefix}/last_interaction_in_nu"].to_numpy()[
                    emitted_packet_mask
                ],
                "Hz",
            ),
            "last_interaction_in_r": u.Quantity(
                hdf[f"{packet_prefix}/last_interaction_in_r"].to_numpy()[
                    emitted_packet_mask
                ],
                "cm",
            ),
            "packet_nus": u.Quantity(
                hdf[f"{packet_prefix}/output_nu"].to_numpy()[
                    emitted_packet_mask
                ],
                "Hz",
            ),
            "packet_energies": u.Quantity(
                hdf[f"{packet_prefix}/output_energy"].to_numpy()[
                    emitted_packet_mask
                ],
                "erg",
            ),
        }
    packet_data["packets_df"] = pd.DataFrame(packet_data)
    process_line_interactions(packet_data, lines_df)
    return packet_data


def parse_species_list_util(species_list):
    """
    Parse user-requested species list and create list of species IDs to be used.

    The function interprets element or ion names and ion ranges, converting them
    into (Z, ion) tuples, where Z is the atomic number and `ion` is the zero-based ionization stage.

    Parameters
    ----------
    species_list : list of str
        List of species that the user wants to show in distinct colors.
        Species can be given as:
        - An ion (e.g. 'Fe II')
        - An element (e.g. 'Ca')
        - A range of ions (e.g. 'Si I - V')
        - A combination of the above (e.g. ['Si II', 'Fe I - III', 'Ca'])

    Returns
    -------
    species_mapped_result : dict
        Dictionary mapping (Z, ion) to lists of (Z, ion) tuples.
    species_list_result : list of tuple
        Flattened list of all (Z, ion) tuples to be used.
    keep_colour_result : list of int
        Atomic numbers of elements that should be grouped by color.
    full_species_list : list of str
        Expanded list of user-requested species in string format.

    Examples
    --------
    'Fe II'        -> [(26, 1)]
    'Ca'           -> [(20, 0), (20, 1), ..., (20, 19)]
    'Si I - V'     -> [(14, 0), (14, 1), (14, 2), (14, 3), (14, 4)]
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
            species_mapped = {}
            for species in species_list:
                # check if a hyphen is present. If it is, then it indicates a
                # range of ions. Add each ion in that range to the list as a new entry
                if "-" in species:
                    # split the string on spaces. First thing in the list is then the element
                    parts = species.split(" ")
                    element = parts[0]
                    ion_range = parts[-1]
                    # Next thing is the ion range
                    # convert the requested ions into numerals
                    range_parts = [
                        part.strip() for part in ion_range.split("-")
                    ]
                    first_ion_numeral = roman_to_int(range_parts[0])
                    second_ion_numeral = roman_to_int(range_parts[-1])
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
            requested_species_ids = []
            keep_colour = []

            # go through each of the requested species. Check whether it is
            # an element or ion (ions have spaces). If it is an element,
            # add all possible ions to the ions list. Otherwise just add
            # the requested ion
            for species in full_species_list:
                if " " in species:
                    species_id = (
                        species_string_to_tuple(species)[0],
                        species_string_to_tuple(species)[1],
                    )
                    requested_species_ids.append([species_id])
                    species_mapped[species_id] = [species_id]
                else:
                    atomic_number = element_symbol2atomic_number(species)
                    species_ids = [
                        (atomic_number, ion_number)
                        for ion_number in np.arange(atomic_number)
                    ]
                    requested_species_ids.append(species_ids)
                    species_mapped[(atomic_number, 0)] = species_ids
                    # add the atomic number to a list so you know that this element should
                    # have all species in the same colour, i.e. it was requested like
                    # species_list = [Si]
                    keep_colour.append(atomic_number)
            requested_species_ids = [
                species_id
                for temp_list in requested_species_ids
                for species_id in temp_list
            ]
            species_mapped_result = species_mapped
            species_list_result = requested_species_ids
            keep_colour_result = keep_colour
    else:
        species_list_result = None
        species_mapped_result = None
        keep_colour_result = None
        full_species_list = None

    return (
        species_mapped_result,
        species_list_result,
        keep_colour_result,
        full_species_list,
    )


def make_colorbar_labels(species, species_list=None, species_mapped=None):
    """
    Generate labels for the colorbar based on species.

    Parameters
    ----------
    species : list of int
        List of species identifiers (Z * 100 + ion) or atomic numbers.
    species_list : list, optional
        Optional list of species to filter against.
    species_mapped : dict, optional
        Mapping from species key (Z * 100 + ion) to lists of species IDs.

    Returns
    -------
    list of str
        List of formatted species labels
    """
    if species_list is None:
        species_name = [
            atomic_number2element_symbol(atomic_num) for atomic_num in species
        ]
    else:
        species_name = []
        for species_key, species_ids in species_mapped.items():
            if any(spec_id in species for spec_id in species_ids):
                if species_key % 100 == 0:
                    label = atomic_number2element_symbol(species_key // 100)
                else:
                    atomic_number = species_key // 100
                    ion_number = species_key % 100
                    ion_numeral = int_to_roman(ion_number + 1)
                    label = f"{atomic_number2element_symbol(atomic_number)} {ion_numeral}"
                species_name.append(label)

    return species_name


def get_spectrum_data(packets_mode, sim):
    """
    Get spectrum data from simulation based on mode.

    Parameters
    ----------
    packets_mode : str
        Packet mode to extract spectrum data for (e.g., 'real', 'virtual').
    sim : Simulation
        Simulation object containing the spectrum solver and data.

    Returns
    -------
    dict
        Dictionary containing:
        - "spectrum_delta_frequency" : Quantity
            Frequency bin width in Hz.
        - "spectrum_frequency_bins" : Quantity
            Frequency bin edges in Hz.
        - "spectrum_luminosity_density_lambda" : Quantity
            Luminosity density in erg / s / Å.
        - "spectrum_wavelength" : Quantity
            Wavelength values in Å.
    """
    packets_type = f"spectrum_{packets_mode}_packets"

    return {
        "spectrum_delta_frequency": getattr(
            sim.spectrum_solver, packets_type
        ).delta_frequency,
        "spectrum_frequency_bins": getattr(
            sim.spectrum_solver, packets_type
        )._frequency,
        "spectrum_luminosity_density_lambda": getattr(
            sim.spectrum_solver, packets_type
        ).luminosity_density_lambda,
        "spectrum_wavelength": getattr(
            sim.spectrum_solver, packets_type
        ).wavelength,
    }


def extract_spectrum_data_hdf(hdf, packets_mode):
    """
    Extract spectrum data from HDF5.

    Parameters
    ----------
    hdf : h5py.File
        Open HDF5 file containing simulation output.
    packets_mode : str
        Packet mode to extract spectrum data for (e.g., 'real', 'virtual').

    Returns
    -------
    dict
        Dictionary containing:
        - "spectrum_delta_frequency" : Quantity
            Frequency bin width in Hz.
        - "spectrum_frequency_bins" : Quantity
            Frequency bin edges in Hz.
        - "spectrum_luminosity_density_lambda" : Quantity
            Luminosity density in erg / s / Å.
        - "spectrum_wavelength" : Quantity
            Wavelength values in Å.
    """
    spectrum_prefix = (
        f"/simulation/spectrum_solver/spectrum_{packets_mode}_packets"
    )
    return {
        "spectrum_delta_frequency": u.Quantity(
            hdf[f"{spectrum_prefix}/scalars"].delta_frequency, "Hz"
        ),
        "spectrum_frequency_bins": u.Quantity(
            hdf[f"{spectrum_prefix}/_frequency"].to_numpy(), "Hz"
        ),
        "spectrum_luminosity_density_lambda": u.Quantity(
            hdf[f"{spectrum_prefix}/luminosity_density_lambda"].to_numpy(),
            "erg / s cm",
        ).to("erg / s AA"),
        "spectrum_wavelength": u.Quantity(
            hdf[f"{spectrum_prefix}/wavelength"].to_numpy(), "cm"
        ).to("AA"),
    }


def get_packet_data(transport_state, packets_mode):
    """Get packet data from transport state based on mode."""
    if packets_mode == "virtual":
        vpacket_tracker = transport_state.vpacket_tracker
        return {
            "last_interaction_type": vpacket_tracker.last_interaction_type,
            "last_line_interaction_in_id": vpacket_tracker.last_interaction_in_id,
            "last_line_interaction_out_id": vpacket_tracker.last_interaction_out_id,
            "last_line_interaction_in_nu": vpacket_tracker.last_interaction_in_nu,
            "last_interaction_in_r": vpacket_tracker.last_interaction_in_r,
            "nus": u.Quantity(vpacket_tracker.nus, "Hz"),
            "energies": u.Quantity(vpacket_tracker.energies, "erg"),
            "lambdas": u.Quantity(vpacket_tracker.nus, "Hz").to(
                "angstrom", u.spectral()
            ),
        }
    # real packets
    mask = transport_state.emitted_packet_mask
    packet_nus = u.Quantity(
        transport_state.packet_collection.output_nus[mask], u.Hz
    )
    return {
        "last_interaction_type": transport_state.last_interaction_type[mask],
        "last_line_interaction_in_id": transport_state.last_line_interaction_in_id[
            mask
        ],
        "last_line_interaction_out_id": transport_state.last_line_interaction_out_id[
            mask
        ],
        "last_line_interaction_in_nu": transport_state.last_interaction_in_nu[
            mask
        ],
        "last_interaction_in_r": transport_state.last_interaction_in_r[mask],
        "nus": packet_nus,
        "energies": transport_state.packet_collection.output_energies[mask],
        "lambdas": packet_nus.to("angstrom", u.spectral()),
    }


def process_line_interactions(packet_data, lines_df):
    """Process line interactions and create line interaction dataframe for both packet modes."""
    for packets_mode in ["real", "virtual"]:
        packets_df = packet_data[packets_mode]["packets_df"]

        if packets_df is not None:
            # Create dataframe of packets that experience line interaction
            line_mask = (packets_df["last_interaction_type"] > -1) & (
                packets_df["last_line_interaction_in_id"] > -1
            )
            packet_data[packets_mode]["packets_df_line_interaction"] = (
                packets_df.loc[line_mask].copy()
            )

            # Add columns for atomic number of last interaction out
            packet_data[packets_mode]["packets_df_line_interaction"][
                "last_line_interaction_atom"
            ] = (
                lines_df["atomic_number"]
                .iloc[
                    packet_data[packets_mode]["packets_df_line_interaction"][
                        "last_line_interaction_out_id"
                    ]
                ]
                .to_numpy()
            )

            # Add columns for the species ID of last interaction
            packet_data[packets_mode]["packets_df_line_interaction"][
                "last_line_interaction_species"
            ] = (
                lines_df["atomic_number"]
                .iloc[
                    packet_data[packets_mode]["packets_df_line_interaction"][
                        "last_line_interaction_out_id"
                    ]
                ]
                .to_numpy()
                * 100
                + lines_df["ion_number"]
                .iloc[
                    packet_data[packets_mode]["packets_df_line_interaction"][
                        "last_line_interaction_out_id"
                    ]
                ]
                .to_numpy()
            )


def extract_packet_data_hdf(hdf, packets_mode):
    """Extract packet data from HDF."""
    if packets_mode == "virtual":
        packet_prefix = "/simulation/transport/transport_state/virt_packet"
        return {
            "last_interaction_type": hdf[
                f"{packet_prefix}_last_interaction_type"
            ],
            "last_line_interaction_in_id": hdf[
                f"{packet_prefix}_last_line_interaction_in_id"
            ],
            "last_line_interaction_out_id": hdf[
                f"{packet_prefix}_last_line_interaction_out_id"
            ],
            "last_line_interaction_in_nu": u.Quantity(
                hdf[f"{packet_prefix}_last_interaction_in_nu"].to_numpy(), "Hz"
            ),
            "last_interaction_in_r": u.Quantity(
                hdf[f"{packet_prefix}_last_interaction_in_r"].to_numpy(), "cm"
            ),
            "packet_nus": u.Quantity(
                hdf[f"{packet_prefix}_nus"].to_numpy(), "Hz"
            ),
            "packet_energies": u.Quantity(
                hdf[f"{packet_prefix}_energies"].to_numpy(), "erg"
            ),
        }
    else:  # real packets
        emitted_packet_mask = hdf[
            "/simulation/transport/transport_state/emitted_packet_mask"
        ].to_numpy()
        packet_prefix = "/simulation/transport/transport_state"
        return {
            "last_interaction_type": hdf[
                f"{packet_prefix}/last_interaction_type"
            ].to_numpy()[emitted_packet_mask],
            "last_line_interaction_in_id": hdf[
                f"{packet_prefix}/last_line_interaction_in_id"
            ].to_numpy()[emitted_packet_mask],
            "last_line_interaction_out_id": hdf[
                f"{packet_prefix}/last_line_interaction_out_id"
            ].to_numpy()[emitted_packet_mask],
            "last_line_interaction_in_nu": u.Quantity(
                hdf[f"{packet_prefix}/last_interaction_in_nu"].to_numpy()[
                    emitted_packet_mask
                ],
                "Hz",
            ),
            "last_interaction_in_r": u.Quantity(
                hdf[f"{packet_prefix}/last_interaction_in_r"].to_numpy()[
                    emitted_packet_mask
                ],
                "cm",
            ),
            "packet_nus": u.Quantity(
                hdf[f"{packet_prefix}/output_nu"].to_numpy()[
                    emitted_packet_mask
                ],
                "Hz",
            ),
            "packet_energies": u.Quantity(
                hdf[f"{packet_prefix}/output_energy"].to_numpy()[
                    emitted_packet_mask
                ],
                "erg",
            ),
        }


def parse_species_list_util(species_list):
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
            species_mapped = {}
            for species in species_list:
                # check if a hyphen is present. If it is, then it indicates a
                # range of ions. Add each ion in that range to the list as a new entry
                if "-" in species:
                    # split the string on spaces. First thing in the list is then the element
                    parts = species.split(" ")
                    element = parts[0]
                    ion_range = parts[-1]
                    # Next thing is the ion range
                    # convert the requested ions into numerals
                    range_parts = [
                        part.strip() for part in ion_range.split("-")
                    ]
                    first_ion_numeral = roman_to_int(range_parts[0])
                    second_ion_numeral = roman_to_int(range_parts[-1])
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
            requested_species_ids = []
            keep_colour = []

            # go through each of the requested species. Check whether it is
            # an element or ion (ions have spaces). If it is an element,
            # add all possible ions to the ions list. Otherwise just add
            # the requested ion
            for species in full_species_list:
                if " " in species:
                    species_id = (
                        species_string_to_tuple(species)[0],
                        species_string_to_tuple(species)[1],
                    )
                    requested_species_ids.append([species_id])
                    species_mapped[species_id] = [species_id]
                else:
                    atomic_number = element_symbol2atomic_number(species)
                    species_ids = [
                        (atomic_number, ion_number)
                        for ion_number in np.arange(atomic_number)
                    ]
                    requested_species_ids.append(species_ids)
                    species_mapped[(atomic_number, 0)] = species_ids
                    # add the atomic number to a list so you know that this element should
                    # have all species in the same colour, i.e. it was requested like
                    # species_list = [Si]
                    keep_colour.append(atomic_number)
            requested_species_ids = [
                species_id
                for temp_list in requested_species_ids
                for species_id in temp_list
            ]
            species_mapped_result = species_mapped
            species_list_result = requested_species_ids
            keep_colour_result = keep_colour
    else:
        species_list_result = None

    return (
        species_mapped_result,
        species_list_result,
        keep_colour_result,
        full_species_list,
    )


def make_colorbar_labels(species, species_list=None, species_mapped=None):
    """
    Generate labels for the colorbar based on species.

    If a species list is provided, uses that to generate labels.
    Otherwise, generates labels from the species in the model.
    """
    if species_list is None:
        species_name = [
            atomic_number2element_symbol(atomic_num) for atomic_num in species
        ]
    else:
        species_name = []
        for species_key, species_ids in species_mapped.items():
            if any(spec_id in species for spec_id in species_ids):
                if species_key % 100 == 0:
                    label = atomic_number2element_symbol(species_key // 100)
                else:
                    atomic_number = species_key // 100
                    ion_number = species_key % 100
                    ion_numeral = int_to_roman(ion_number + 1)
                    label = f"{atomic_number2element_symbol(atomic_number)} {ion_numeral}"
                species_name.append(label)

    return species_name


def get_spectrum_data(packets_mode, sim):
    """Get spectrum data from simulation based on mode."""
    packets_type = f"spectrum_{packets_mode}_packets"

    return {
        "spectrum_delta_frequency": getattr(
            sim.spectrum_solver, packets_type
        ).delta_frequency,
        "spectrum_frequency_bins": getattr(
            sim.spectrum_solver, packets_type
        )._frequency,
        "spectrum_luminosity_density_lambda": getattr(
            sim.spectrum_solver, packets_type
        ).luminosity_density_lambda,
        "spectrum_wavelength": getattr(
            sim.spectrum_solver, packets_type
        ).wavelength,
    }


def extract_spectrum_data_hdf(hdf, packets_mode):
    """Extract spectrum data from HDF."""
    spectrum_prefix = (
        f"/simulation/spectrum_solver/spectrum_{packets_mode}_packets"
    )
    return {
        "spectrum_delta_frequency": u.Quantity(
            hdf[f"{spectrum_prefix}/scalars"].delta_frequency, "Hz"
        ),
        "spectrum_frequency_bins": u.Quantity(
            hdf[f"{spectrum_prefix}/_frequency"].to_numpy(), "Hz"
        ),
        "spectrum_luminosity_density_lambda": u.Quantity(
            hdf[f"{spectrum_prefix}/luminosity_density_lambda"].to_numpy(),
            "erg / s cm",
        ).to("erg / s AA"),
        "spectrum_wavelength": u.Quantity(
            hdf[f"{spectrum_prefix}/wavelength"].to_numpy(), "cm"
        ).to("AA"),
    }
