"""Utility functions to be used in plotting."""

import re

import numpy as np

from tardis.util.base import (
    element_symbol2atomic_number,
    int_to_roman,
    roman_to_int,
    species_string_to_tuple,
)
from tardis.visualization.tools.visualization_data import (
    VisualizationData,
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

def create_packet_data_dict_from_simulation(sim):
    """
    Create a dictionary containing virtual and real packet data based on simulation state.

    Parameters
    ----------
    sim : tardis.simulation.Simulation
        TARDIS Simulation object produced by running a simulation

    Returns
    -------
    dict
        Dictionary containing 'virtual' and 'real' SimulationPacketData instances
    """
    packet_data = {
        "real": VisualizationData.from_simulation(sim, "real")
    }
    if sim.transport.transport_state.virt_logging:
        packet_data["virtual"] = VisualizationData.from_simulation(sim, "virtual")
    else:
        packet_data["virtual"] = None

    return packet_data

def create_packet_data_dict_from_hdf(hdf_fpath, packets_mode=None):
    """
    Create a dictionary containing virtual and real packet data from HDF file.

    Parameters
    ----------
    hdf_fpath : str
        Valid path to the HDF file where simulation is saved
    packets_mode : {'virtual', 'real', None}
        Mode of packets to be considered. If None, both modes are returned.

    Returns
    -------
    dict
        Dictionary containing 'virtual' and 'real' SimulationPacketData instances
    """
    if packets_mode not in [None, "virtual", "real"]:
        raise ValueError(
            "Invalid value passed to packets_mode. Only "
            "allowed values are 'virtual', 'real' or None"
        )
    if packets_mode == "virtual":
        return {
            "virtual": VisualizationData.from_hdf(hdf_fpath, "virtual"),
            "real": None
        }
    if packets_mode == "real":
        return {
            "virtual": None,
            "real": VisualizationData.from_hdf(hdf_fpath, "real")
        }
    return {
        "virtual": VisualizationData.from_hdf(hdf_fpath, "virtual"),
        "real": VisualizationData.from_hdf(hdf_fpath, "real")
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

    Returns
    -------
    dict
        A dictionary containing:
        - full_species_list: List of expanded species (e.g. Si I - V -> [Si I, Si II, ...]).
        - species_mapped: Mapping of species ids to species names.
        - keep_colour: List of atomic numbers to group elements with consistent colors.
    """
    if species_list is None:
        return {
        "full_species_list": None,
        "species_mapped": None,
        "keep_colour": None,
        "species_list": None,
    }


    if any(char.isdigit() for char in " ".join(species_list)):
        raise ValueError("All species must be in Roman numeral form, e.g., Si II")

    full_species_list = []
    species_mapped = {}
    keep_colour = []

    for species in species_list:
        if "-" in species:
            element = species.split(" ")[0]
            first_ion_numeral = roman_to_int(species.split(" ")[-1].split("-")[0])
            second_ion_numeral = roman_to_int(species.split(" ")[-1].split("-")[-1])

            for ion_number in range(first_ion_numeral, second_ion_numeral + 1):
                full_species_list.append(f"{element} {int_to_roman(ion_number)}")
        else:
            full_species_list.append(species)

    requested_species_ids = []

    for species in full_species_list:
        if " " in species:
            species_id = (
                species_string_to_tuple(species)[0] * 100
                + species_string_to_tuple(species)[1]
            )
            requested_species_ids.append([species_id])
            species_mapped[species_id] = [species_id]
        else:
            atomic_number = element_symbol2atomic_number(species)
            species_ids = [
                atomic_number * 100 + ion_number for ion_number in range(atomic_number)
            ]
            requested_species_ids.append(species_ids)
            species_mapped[atomic_number * 100] = species_ids
            keep_colour.append(atomic_number)

    requested_species_ids = [
        species_id for temp_list in requested_species_ids for species_id in temp_list
    ]

    return {
        "full_species_list": full_species_list,
        "species_mapped": species_mapped,
        "keep_colour": keep_colour,
        "species_list": requested_species_ids,
    }
