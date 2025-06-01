import astropy.units as u
import numpy as np
import pandas as pd


def calculate_interaction_luminosities(
    transport_state,
    spectrum_solver,
    plasma,
    packets_mode,
    packet_wvl_range,
    species_list=None,
    distance=None,
    mode="absorption",
):
    """
    Calculate line interaction luminosities.

    Parameters
    ----------
    sim : object
        Simulation instance with transport data.
    packets_mode : str
        'virtual' or 'real' packet selection mode.
    packet_wvl_range : astropy.units.Quantity or None
        Wavelength range mask for packets.
    species_list : array-like or None
        Species identifiers for grouping luminosities.
    distance : astropy.units.Quantity or None
        Distance for flux scaling. If None, no flux scaling is done.
    mode : {'absorption', 'emission'}
        Mode of interaction to compute.

    Returns
    -------
    luminosities_df : pandas.DataFrame
        Luminosity values binned over wavelength.
    species_present : numpy.ndarray
        Identifiers of species present in the data.
    """
    if packets_mode == "virtual":
        packets_df = pd.DataFrame(
            {
                "nus": transport_state.vpacket_tracker.nus,
                "energies": transport_state.vpacket_tracker.energies,
                "last_interaction_type": transport_state.vpacket_tracker.last_interaction_type,
                "last_line_interaction_in_id": transport_state.vpacket_tracker.last_interaction_in_id,
                "last_line_interaction_in_nu": transport_state.vpacket_tracker.last_interaction_in_nu,
                "last_line_interaction_out_id": transport_state.vpacket_tracker.last_interaction_out_id,
            }
        )
        spectrum_packets = spectrum_solver.spectrum_virtual_packets
    else:
        emitted_mask = transport_state.emitted_packet_mask
        packets_df = pd.DataFrame(
            {
                "nus": transport_state.packet_collection.output_nus[emitted_mask],
                "energies": transport_state.packet_collection.output_energies[
                    emitted_mask
                ],
                "last_interaction_type": transport_state.last_interaction_type[
                    emitted_mask
                ],
                "last_line_interaction_in_id": transport_state.last_line_interaction_in_id[
                    emitted_mask
                ],
                "last_line_interaction_in_nu": transport_state.last_interaction_in_nu[
                    emitted_mask
                ],
                "last_line_interaction_out_id": transport_state.last_line_interaction_out_id[
                    emitted_mask
                ],
            }
        )
        spectrum_packets = spectrum_solver.spectrum_real_packets

    spectrum_delta_frequency = spectrum_packets.delta_frequency
    plot_wavelength = spectrum_packets.wavelength
    plot_frequency = spectrum_packets._frequency[:-1]
    plot_frequency_bins = spectrum_packets._frequency
    time_of_simulation = transport_state.packet_collection.time_of_simulation
    lines_df = plasma.atomic_data.lines.reset_index().set_index("line_id")

    line_mask = (packets_df["last_interaction_type"] > -1) & (
        packets_df["last_line_interaction_in_id"] > -1
    )
    packets_df_line_interaction = packets_df.loc[line_mask].copy()
    packets_df_line_interaction["last_line_interaction_atom"] = (
        lines_df["atomic_number"]
        .iloc[packets_df_line_interaction["last_line_interaction_out_id"]]
        .to_numpy()
    )
    packets_df_line_interaction["last_line_interaction_species"] = (
        lines_df["atomic_number"]
        .iloc[packets_df_line_interaction["last_line_interaction_out_id"]]
        .to_numpy()
        * 100
        + lines_df["ion_number"]
        .iloc[packets_df_line_interaction["last_line_interaction_out_id"]]
        .to_numpy()
    )

    if mode == "absorption":
        if packet_wvl_range is None:
            packet_nu_line_range_mask = np.ones(
                packets_df_line_interaction.shape[0], dtype=bool
            )
        else:
            packet_nu_range = packet_wvl_range.to("Hz", u.spectral())
            packet_nu_line_range_mask = (
                packets_df_line_interaction["last_line_interaction_in_nu"]
                < packet_nu_range[0]
            ) & (
                packets_df_line_interaction["last_line_interaction_in_nu"]
                > packet_nu_range[1]
            )
    else:
        if packet_wvl_range is None:
            packet_nu_range_mask = np.ones(packets_df.shape[0], dtype=bool)
            packet_nu_line_range_mask = np.ones(
                packets_df_line_interaction.shape[0], dtype=bool
            )
        else:
            packet_nu_range = packet_wvl_range.to("Hz", u.spectral())
            packet_nu_range_mask = (packets_df["nus"] < packet_nu_range[0]) & (
                packets_df["nus"] > packet_nu_range[1]
            )
            packet_nu_line_range_mask = (
                packets_df_line_interaction["nus"] < packet_nu_range[0]
            ) & (packets_df_line_interaction["nus"] > packet_nu_range[1])

    luminosities_df = pd.DataFrame(index=plot_wavelength)
    if distance is None:
        lum_to_flux = 1
    elif distance <= 0:
        raise ValueError(
            "distance passed must be greater than 0. If you intended "
            "to plot luminosities instead of flux, set distance=None "
            "or don't specify distance parameter in the function call."
        )
    else:
        lum_to_flux = 4.0 * np.pi * (distance.to("cm")) ** 2

    if mode == "absorption":
        if species_list is None:
            packets_df_grouped = packets_df_line_interaction.loc[
                packet_nu_line_range_mask
            ].groupby(by="last_line_interaction_atom")
        else:
            packets_df_grouped = packets_df_line_interaction.loc[
                packet_nu_line_range_mask
            ].groupby(by="last_line_interaction_species")
        for identifier, group in packets_df_grouped:
            hist_el = np.histogram(
                group["last_line_interaction_in_nu"],
                bins=plot_frequency_bins.value,
                weights=group["energies"] / lum_to_flux / time_of_simulation,
            )
            L_nu_el = hist_el[0] * u.erg / u.s / spectrum_delta_frequency
            L_lambda_el = L_nu_el * plot_frequency / plot_wavelength
            luminosities_df[identifier] = L_lambda_el.value
        species_present = np.array(list(packets_df_grouped.groups.keys()))
        return luminosities_df, species_present

    weights = (
        packets_df["energies"][packet_nu_range_mask] / lum_to_flux
    ) / time_of_simulation
    mask_noint = packets_df["last_interaction_type"][packet_nu_range_mask] == -1
    hist_noint = np.histogram(
        packets_df["nus"][packet_nu_range_mask][mask_noint],
        bins=plot_frequency_bins.value,
        weights=weights[mask_noint],
        density=False,
    )
    L_nu_noint = hist_noint[0] * u.erg / u.s / spectrum_delta_frequency
    L_lambda_noint = L_nu_noint * plot_frequency / plot_wavelength
    luminosities_df["noint"] = L_lambda_noint.value

    mask_escatter = (packets_df["last_interaction_type"][packet_nu_range_mask] == 1) & (
        packets_df["last_line_interaction_in_id"][packet_nu_range_mask] == -1
    )
    hist_escatter = np.histogram(
        packets_df["nus"][packet_nu_range_mask][mask_escatter],
        bins=plot_frequency_bins.value,
        weights=weights[mask_escatter],
        density=False,
    )
    L_nu_escatter = hist_escatter[0] * u.erg / u.s / spectrum_delta_frequency
    L_lambda_escatter = L_nu_escatter * plot_frequency / plot_wavelength
    luminosities_df["escatter"] = L_lambda_escatter.value

    if species_list is None:
        packets_df_grouped = packets_df_line_interaction.loc[
            packet_nu_line_range_mask
        ].groupby(by="last_line_interaction_atom")
    else:
        packets_df_grouped = packets_df_line_interaction.loc[
            packet_nu_line_range_mask
        ].groupby(by="last_line_interaction_species")
    for identifier, group in packets_df_grouped:
        hist_el = np.histogram(
            group["nus"],
            bins=plot_frequency_bins.value,
            weights=group["energies"] / lum_to_flux / time_of_simulation,
        )
        L_nu_el = hist_el[0] * u.erg / u.s / spectrum_delta_frequency
        L_lambda_el = L_nu_el * plot_frequency / plot_wavelength
        luminosities_df[identifier] = L_lambda_el.value
    species_present = np.array(list(packets_df_grouped.groups.keys()))
    return luminosities_df, species_present


def calculate_absorption_luminosities(
    sim, packets_mode, packet_wvl_range, species_list=None, distance=None
):
    """
    Thin wrapper for absorption luminosities.
    """
    transport_state = sim.transport.transport_state
    spectrum_solver = sim.spectrum_solver
    plasma = sim.plasma
    return calculate_interaction_luminosities(
        transport_state,
        spectrum_solver,
        plasma,
        packets_mode,
        packet_wvl_range,
        species_list,
        distance,
        mode="absorption",
    )


def calculate_emission_luminosities(
    sim, packets_mode, packet_wvl_range, species_list=None, distance=None
):
    """
    Thin wrapper for emission luminosities.
    """
    transport_state = sim.transport.transport_state
    spectrum_solver = sim.spectrum_solver
    plasma = sim.plasma

    return calculate_interaction_luminosities(
        transport_state,
        spectrum_solver,
        plasma,
        packets_mode,
        packet_wvl_range,
        species_list,
        distance,
        mode="emission",
    )
