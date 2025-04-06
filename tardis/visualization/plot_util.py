"""Utility functions to be used in plotting."""

import re

import astropy.units as u
import numpy as np

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


def get_packet_data(transport_state, packets_mode):
    """Get packet data from transport state based on mode."""
    if packets_mode == "virtual":
        vpacket_tracker = transport_state.vpacket_tracker
        return {
            'last_interaction_type': vpacket_tracker.last_interaction_type,
            'last_line_interaction_in_id': vpacket_tracker.last_interaction_in_id,
            'last_line_interaction_out_id': vpacket_tracker.last_interaction_out_id,
            'last_line_interaction_in_nu': vpacket_tracker.last_interaction_in_nu,
            'last_interaction_in_r': vpacket_tracker.last_interaction_in_r,
            'nus': u.Quantity(vpacket_tracker.nus, "Hz"),
            'energies': u.Quantity(vpacket_tracker.energies, "erg"),
            'lambdas': u.Quantity(vpacket_tracker.nus, "Hz").to("angstrom", u.spectral()),
        }
    # real packets
    mask = transport_state.emitted_packet_mask
    packet_nus = u.Quantity(transport_state.packet_collection.output_nus[mask], u.Hz)
    return {
        'last_interaction_type': transport_state.last_interaction_type[mask],
        'last_line_interaction_in_id': transport_state.last_line_interaction_in_id[mask],
        'last_line_interaction_out_id': transport_state.last_line_interaction_out_id[mask],
        'last_line_interaction_in_nu': transport_state.last_interaction_in_nu[mask],
        'last_interaction_in_r': transport_state.last_interaction_in_r[mask],
        'nus': packet_nus,
        'energies': transport_state.packet_collection.output_energies[mask],
        'lambdas': packet_nus.to("angstrom", u.spectral()),
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
            packet_data[packets_mode]["packets_df_line_interaction"] = packets_df.loc[line_mask].copy()

            # Add columns for atomic number of last interaction out
            packet_data[packets_mode]["packets_df_line_interaction"]["last_line_interaction_atom"] = (
                lines_df["atomic_number"]
                .iloc[packet_data[packets_mode]["packets_df_line_interaction"]["last_line_interaction_out_id"]]
                .to_numpy()
            )

            # Add columns for the species ID of last interaction
            packet_data[packets_mode]["packets_df_line_interaction"]["last_line_interaction_species"] = (
                lines_df["atomic_number"]
                .iloc[packet_data[packets_mode]["packets_df_line_interaction"]["last_line_interaction_out_id"]]
                .to_numpy()
                * 100
                + lines_df["ion_number"]
                .iloc[packet_data[packets_mode]["packets_df_line_interaction"]["last_line_interaction_out_id"]]
                .to_numpy()
            )

def extract_packet_data_hdf(hdf, packets_mode):
    """Extract packet data from HDF."""
    if packets_mode == "virtual":
        packet_prefix = "/simulation/transport/transport_state/virt_packet"
        return {
            "last_interaction_type": hdf[f"{packet_prefix}_last_interaction_type"],
            "last_line_interaction_in_id": hdf[f"{packet_prefix}_last_line_interaction_in_id"],
            "last_line_interaction_out_id": hdf[f"{packet_prefix}_last_line_interaction_out_id"],
            "last_line_interaction_in_nu": u.Quantity(hdf[f"{packet_prefix}_last_interaction_in_nu"].to_numpy(), "Hz"),
            "last_interaction_in_r": u.Quantity(hdf[f"{packet_prefix}_last_interaction_in_r"].to_numpy(), "cm"),
            "packet_nus": u.Quantity(hdf[f"{packet_prefix}_nus"].to_numpy(), "Hz"),
            "packet_energies": u.Quantity(hdf[f"{packet_prefix}_energies"].to_numpy(), "erg"),
        }
    else:  # real packets
        emitted_packet_mask = hdf["/simulation/transport/transport_state/emitted_packet_mask"].to_numpy()
        packet_prefix = "/simulation/transport/transport_state"
        return {
            "last_interaction_type": hdf[f"{packet_prefix}/last_interaction_type"].to_numpy()[emitted_packet_mask],
            "last_line_interaction_in_id": hdf[f"{packet_prefix}/last_line_interaction_in_id"].to_numpy()[emitted_packet_mask],
            "last_line_interaction_out_id": hdf[f"{packet_prefix}/last_line_interaction_out_id"].to_numpy()[emitted_packet_mask],
            "last_line_interaction_in_nu": u.Quantity(
                hdf[f"{packet_prefix}/last_interaction_in_nu"].to_numpy()[emitted_packet_mask], "Hz"
            ),
            "last_interaction_in_r": u.Quantity(
                hdf[f"{packet_prefix}/last_interaction_in_r"].to_numpy()[emitted_packet_mask], "cm"
            ),
            "packet_nus": u.Quantity(hdf[f"{packet_prefix}/output_nu"].to_numpy()[emitted_packet_mask], "Hz"),
            "packet_energies": u.Quantity(hdf[f"{packet_prefix}/output_energy"].to_numpy()[emitted_packet_mask], "erg"),
        }