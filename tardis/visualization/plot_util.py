"""Utility functions to be used in plotting."""

import re

import astropy.units as u
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tardis.transport.montecarlo.packets.radiative_packet import InteractionType

from tardis.util.base import (
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


def get_hex_color_strings(length, name="jet"):
    """
    Generate a list of hex color strings from a discrete colormap.

    Parameters
    ----------
    length : int
        Number of discrete colors to extract from the colormap.
    name : str, optional
        Name of the Matplotlib colormap to use (default is 'jet').

    Returns
    -------
    list of str
        List of hex color strings (['#ff0000', '#00ff00', '#0000ff'])
    """
    cmap = plt.get_cmap(name, length)
    return [mcolors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N)]


def extract_and_process_packet_data(simulation, packets_mode, include_shell_id=False):
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

    include_shell_id : bool, optional
        Whether to include shell_id information in the output (default: False)

    Returns
    -------
    dict
        Dictionary containing raw packet data, the full DataFrame `packets_df`,
        and a filtered `packets_df_line_interaction` with line interaction info.
    """
    if hasattr(simulation, "transport_state"): # for workflows
        transport_state = simulation.transport_state
        lines = simulation.plasma_solver.atomic_data.lines
    else:
        transport_state = simulation.transport.transport_state
        lines = simulation.plasma.atomic_data.lines

    lines_df = lines.reset_index().set_index("line_id")

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
        # Get emitted packets that had line interactions
        df = transport_state.tracker_full_df
        
        # Get final interaction per packet where status is EMITTED
        emitted_final = df.groupby(level='packet_id').last()
        emitted_packets = emitted_final[emitted_final['status'] == 'EMITTED'].index
        
        # Get last line interaction for each emitted packet
        line_interactions = df[df['interaction_type'] == 'LINE']
        emitted_line_interactions = line_interactions[
            line_interactions.index.get_level_values('packet_id').isin(emitted_packets)
        ]
        result_df = emitted_line_interactions.groupby(level='packet_id').last()
        
        # Extract packet collection data for these packets
        packet_indices = result_df.index.values
        packet_nus = u.Quantity(
            transport_state.packet_collection.output_nus[packet_indices], u.Hz
        )
        
        packet_data = {
            "last_interaction_type": result_df["interaction_type"].values,
            "last_line_interaction_in_id": result_df["line_absorb_id"].values,
            "last_line_interaction_out_id": result_df["line_emit_id"].values,
            "last_line_interaction_in_nu": result_df["before_nu"].values,
            "last_interaction_in_r": result_df["radius"].values,
            "nus": packet_nus,
            "energies": transport_state.packet_collection.output_energies[packet_indices],
            "lambdas": packet_nus.to("angstrom", u.spectral()),
        }

        if include_shell_id:
            packet_data["last_line_interaction_shell_id"] = result_df["after_shell_id"].values

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
        line_mask = (packets_df["last_interaction_type"] != "NO_INTERACTION") & (
            packets_df["last_line_interaction_in_id"] > -1
        )
        packet_data["packets_df_line_interaction"] = packets_df.loc[
            line_mask
        ].copy()

        # Add columns for atomic number of last interaction out
        packet_data["packets_df_line_interaction"][
            "last_line_interaction_atom"
        ] = list(
            zip(
                lines_df["atomic_number"]
                .iloc[
                    packet_data["packets_df_line_interaction"][
                        "last_line_interaction_out_id"
                    ]
                ]
                .to_numpy(),
                [0] * len(packet_data["packets_df_line_interaction"]),
            )
        )

        # Add columns for the species ID of last interaction
        packet_data["packets_df_line_interaction"][
            "last_line_interaction_species"
        ] = list(zip(
            lines_df["atomic_number"]
            .iloc[
                packet_data["packets_df_line_interaction"][
                    "last_line_interaction_out_id"
                ]
            ]
            .to_numpy(),
            lines_df["ion_number"]
            .iloc[
                packet_data["packets_df_line_interaction"][
                    "last_line_interaction_out_id"
                ]
            ]
            .to_numpy()
        ))


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


def expand_species_list(species_list):
    """
    Expand a species list into a fully-resolved list of species strings.

    This includes:
    - Expanding ion ranges like 'Si I - V' into ['Si I', 'Si II', ..., 'Si V']
    - Keeping individual ions or elements as-is

    Parameters
    ----------
    species_list : list of str
        List of species requested by the user.

    Returns
    -------
    full_species_list : list of str
        Expanded list of species strings in 'Element Ion' format.

    Raises
    ------
    ValueError
        If any digit is found in the input (species must use Roman numerals).
    """
    # check if there are any digits in the species list. If there are, then exit.
    # species_list should only contain species in the Roman numeral
    # format, e.g. Si II, and each ion must contain a space
    if any(char.isdigit() for char in " ".join(species_list)):
        raise ValueError(
            "All species must be in Roman numeral form, e.g. Si II"
        )

    full_species_list = []
    for species in species_list:
        if "-" in species:
            element, ion_range = species.split(" ")
            first_ion_numeral, last_ion_numeral = map(
                roman_to_int, ion_range.partition("-")[::2]
            )
            for ion_number in range(first_ion_numeral, last_ion_numeral + 1):
                full_species_list.append(
                    f"{element} {int_to_roman(ion_number)}"
                )
        else:
            full_species_list.append(species)

    return full_species_list


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
    elements_with_shared_color : list of int
        Atomic numbers of elements that should be grouped by color.
    full_species_list : list of str
        Expanded list of user-requested species in string format.

    Examples
    --------
    'Fe II'        -> [(26, 1)]
    'Ca'           -> [(20, 0), (20, 1), ..., (20, 19)]
    'Si I-V'       -> [(14, 0), (14, 1), (14, 2), (14, 3), (14, 4)]
    """
    if species_list is None:
        return None, None, None, None

    full_species_list = expand_species_list(species_list)
    requested_species_ids = []
    elements_with_shared_color = []
    species_mapped = {}

    # go through each of the requested species. Check whether it is
    # an element or ion (ions have spaces). If it is an element,
    # add all possible ions to the ions list. Otherwise just add
    # the requested ion
    for species in full_species_list:
        if " " in species:
            species_id = species_string_to_tuple(species)
            requested_species_ids.append([species_id])
            species_mapped[species_id] = [species_id]
        else:
            atomic_number = element_symbol2atomic_number(species)
            species_ids = [
                (atomic_number, ion_number)
                for ion_number in range(atomic_number)
            ]
            requested_species_ids.append(species_ids)
            species_mapped[(atomic_number, 0)] = species_ids
            # add the atomic number to a list so you know that this element should
            # have all species in the same colour, i.e. it was requested like
            # species_list = [Si]
            elements_with_shared_color.append(atomic_number)

    species_list_result = [
        species_id for group in requested_species_ids for species_id in group
    ]

    return (
        species_mapped,
        species_list_result,
        elements_with_shared_color,
        full_species_list,
    )


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


def create_wavelength_mask(
    packet_data, packets_mode, packet_wvl_range, df_key, column_name
):
    """
    Create mask for packets based on wavelength range.

    Parameters
    ----------
    packet_data : dict or pd.DataFrame
        Either a nested dict with packet data or a DataFrame directly
    packets_mode : str or None
        'virtual' or 'real' packets mode. Required if packet_data is a dict, ignored if DataFrame
    packet_wvl_range : astropy.Quantity or None
        Wavelength range to filter packets
    df_key : str or None
        Key for the dataframe in packet_data ('packets_df' or 'packets_df_line_interaction').
        Required if packet_data is a dict, ignored if DataFrame
    column_name : str
        Column name to filter on ('nus' or 'last_line_interaction_in_nu')

    Returns
    -------
    np.array
        Boolean mask for packets in the specified wavelength range
    """
    if isinstance(packet_data, pd.DataFrame):
        df = packet_data
    else:
        df = packet_data[packets_mode][df_key]

    if packet_wvl_range is None:
        return np.ones(df.shape[0], dtype=bool)

    packet_nu_range = packet_wvl_range.to("Hz", u.spectral())

    return (df[column_name] < packet_nu_range[0]) & (
        df[column_name] > packet_nu_range[1]
    )
