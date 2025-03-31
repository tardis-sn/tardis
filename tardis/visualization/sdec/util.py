import astropy.units as u
import numpy as np
import pandas as pd


def calculate_absorption_luminosities(
    sim, packets_mode, packet_wvl_range, species_list=None, distance=None
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
    transport_state = sim.transport.transport_state

    # Set up appropriate data sources based on packet mode
    if packets_mode == "virtual":
        # Virtual packets data
        packets_df = pd.DataFrame({
            "nus": transport_state.vpacket_tracker.nus,
            "energies": transport_state.vpacket_tracker.energies,
            "last_interaction_type": transport_state.vpacket_tracker.last_interaction_type,
            "last_line_interaction_in_id": transport_state.vpacket_tracker.last_interaction_in_id,
            "last_line_interaction_in_nu": transport_state.vpacket_tracker.last_interaction_in_nu,
            "last_line_interaction_out_id": transport_state.vpacket_tracker.last_interaction_out_id
        })
        spectrum_delta_frequency = sim.spectrum_solver.spectrum_virtual_packets.delta_frequency
        plot_wavelength = sim.spectrum_solver.spectrum_virtual_packets.wavelength
        plot_frequency = sim.spectrum_solver.spectrum_virtual_packets._frequency[:-1]
        plot_frequency_bins = sim.spectrum_solver.spectrum_virtual_packets._frequency

    else:  # real packets
        emitted_mask = transport_state.emitted_packet_mask
        packets_df = pd.DataFrame({
            "nus": transport_state.packet_collection.output_nus[emitted_mask],
            "energies": transport_state.packet_collection.output_energies[emitted_mask],
            "last_interaction_type": transport_state.last_interaction_type[emitted_mask],
            "last_line_interaction_in_id": transport_state.last_line_interaction_in_id[emitted_mask],
            "last_line_interaction_in_nu": transport_state.last_interaction_in_nu[
                    transport_state.emitted_packet_mask
                ],
            "last_line_interaction_out_id": transport_state.last_line_interaction_out_id[emitted_mask]
        })
        spectrum_delta_frequency = sim.spectrum_solver.spectrum_real_packets.delta_frequency
        plot_wavelength = sim.spectrum_solver.spectrum_real_packets.wavelength
        plot_frequency = sim.spectrum_solver.spectrum_real_packets._frequency[:-1]
        plot_frequency_bins = sim.spectrum_solver.spectrum_real_packets._frequency

    time_of_simulation = transport_state.packet_collection.time_of_simulation
    lines_df = sim.plasma.atomic_data.lines.reset_index().set_index("line_id")

    # Create packets_df_line_interaction
    line_mask = (packets_df["last_interaction_type"] > -1) & (packets_df["last_line_interaction_in_id"] > -1)
    packets_df_line_interaction = packets_df.loc[line_mask].copy()

    # Add atomic number and species columns
    packets_df_line_interaction["last_line_interaction_atom"] = (
        lines_df["atomic_number"]
        .iloc[packets_df_line_interaction["last_line_interaction_out_id"]]
        .to_numpy()
    )
    packets_df_line_interaction["last_line_interaction_species"] = (
        lines_df["atomic_number"]
        .iloc[packets_df_line_interaction["last_line_interaction_out_id"]]
        .to_numpy() * 100
        + lines_df["ion_number"]
        .iloc[packets_df_line_interaction["last_line_interaction_out_id"]]
        .to_numpy()
    )
    if packet_wvl_range is None:
        packet_nu_line_range_mask = np.ones(
            packets_df_line_interaction.shape[0],
            dtype=bool,
        )
    else:
        packet_nu_range = packet_wvl_range.to("Hz", u.spectral())
        packet_nu_line_range_mask = (
            packets_df_line_interaction[
                "last_line_interaction_in_nu"
            ]
            < packet_nu_range[0]
        ) & (
            packets_df_line_interaction[
                "last_line_interaction_in_nu"
            ]
            > packet_nu_range[1]
        )

    luminosities_df = pd.DataFrame(index=plot_wavelength)
         # Calculate the area term to convert luminosity to flux
    if distance is None:
        lum_to_flux = 1  # so that this term will have no effect
    else:
        if distance <= 0:
            raise ValueError(
                "distance passed must be greater than 0. If you intended "
                "to plot luminosities instead of flux, set distance=None "
                "or don't specify distance parameter in the function call."
            )
        else:
            lum_to_flux = 4.0 * np.pi * (distance.to("cm")) ** 2
    # Group packets_df by atomic number of elements with which packets
    # had their last absorption (interaction in)
    # or if species_list is requested then group by species id
    if species_list is None:
        packets_df_grouped = (
            packets_df_line_interaction.loc[packet_nu_line_range_mask]
            .groupby(by="last_line_interaction_atom")
        )
    else:
        packets_df_grouped = (
            packets_df_line_interaction.loc[packet_nu_line_range_mask]
            .groupby(by="last_line_interaction_species")
        )

    for identifier, group in packets_df_grouped:
        # Histogram of specific species
        hist_el = np.histogram(
            group["last_line_interaction_in_nu"],
            bins=plot_frequency_bins.value,
            weights=group["energies"]
            / lum_to_flux
            / time_of_simulation,
        )

        # Convert to luminosity density lambda
        L_nu_el = (
            hist_el[0]
            * u.erg
            / u.s
            / spectrum_delta_frequency
        )
        L_lambda_el = L_nu_el * plot_frequency / plot_wavelength

        luminosities_df[identifier] = L_lambda_el.value

    absorption_species_present = np.array(
        list(packets_df_grouped.groups.keys())
    )

    return luminosities_df, absorption_species_present


def calculate_emission_luminosities(sim, packets_mode, packet_wvl_range, species_list=None, distance=None):
    """
    Calculate luminosities for the emission part of SDEC plot.

    Parameters
    ----------
    sim : tardis.simulation.Simulation
        TARDIS Simulation object produced by running a simulation
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
    transport_state = sim.transport.transport_state

    # Set up appropriate data sources based on packet mode
    if packets_mode == "virtual":
        # Virtual packets data
        packets_df = pd.DataFrame({
            "nus": transport_state.vpacket_tracker.nus,
            "energies": transport_state.vpacket_tracker.energies,
            "last_interaction_type": transport_state.vpacket_tracker.last_interaction_type,
            "last_line_interaction_in_id": transport_state.vpacket_tracker.last_interaction_in_id,
            "last_line_interaction_in_nu": transport_state.vpacket_tracker.last_interaction_in_nu,
            "last_line_interaction_out_id": transport_state.vpacket_tracker.last_interaction_out_id
        })
        spectrum_delta_frequency = sim.spectrum_solver.spectrum_virtual_packets.delta_frequency
        plot_wavelength = sim.spectrum_solver.spectrum_virtual_packets.wavelength
        plot_frequency = sim.spectrum_solver.spectrum_virtual_packets._frequency[:-1]
        plot_frequency_bins = sim.spectrum_solver.spectrum_virtual_packets._frequency

    else:  # real packets
        emitted_mask = transport_state.emitted_packet_mask
        packets_df = pd.DataFrame({
            "nus": transport_state.packet_collection.output_nus[emitted_mask],
            "energies": transport_state.packet_collection.output_energies[emitted_mask],
            "last_interaction_type": transport_state.last_interaction_type[emitted_mask],
            "last_line_interaction_in_id": transport_state.last_line_interaction_in_id[emitted_mask],
            "last_line_interaction_in_nu": transport_state.last_interaction_in_nu[
                    transport_state.emitted_packet_mask
                ],
            "last_line_interaction_out_id": transport_state.last_line_interaction_out_id[emitted_mask]
        })
        spectrum_delta_frequency = sim.spectrum_solver.spectrum_real_packets.delta_frequency
        plot_wavelength = sim.spectrum_solver.spectrum_real_packets.wavelength
        plot_frequency = sim.spectrum_solver.spectrum_real_packets._frequency[:-1]
        plot_frequency_bins = sim.spectrum_solver.spectrum_real_packets._frequency

    time_of_simulation = transport_state.packet_collection.time_of_simulation
    lines_df = sim.plasma.atomic_data.lines.reset_index().set_index("line_id")

    # Create packets_df_line_interaction
    line_mask = (packets_df["last_interaction_type"] > -1) & (packets_df["last_line_interaction_in_id"] > -1)
    packets_df_line_interaction = packets_df.loc[line_mask].copy()

    # Add atomic number and species columns
    packets_df_line_interaction["last_line_interaction_atom"] = (
        lines_df["atomic_number"]
        .iloc[packets_df_line_interaction["last_line_interaction_out_id"]]
        .to_numpy()
    )
    packets_df_line_interaction["last_line_interaction_species"] = (
        lines_df["atomic_number"]
        .iloc[packets_df_line_interaction["last_line_interaction_out_id"]]
        .to_numpy() * 100
        + lines_df["ion_number"]
        .iloc[packets_df_line_interaction["last_line_interaction_out_id"]]
        .to_numpy()
    )

    # Calculate masks to be applied on packets data based on packet_wvl_range
    if packet_wvl_range is None:
        packet_nu_range_mask = np.ones(
            packets_df.shape[0], dtype=bool
        )
        packet_nu_line_range_mask = np.ones(
            packets_df_line_interaction.shape[0],
            dtype=bool,
        )
    else:
        packet_nu_range = packet_wvl_range.to("Hz", u.spectral())
        packet_nu_range_mask = (
            packets_df["nus"] < packet_nu_range[0]
        ) & (packets_df["nus"] > packet_nu_range[1])
        packet_nu_line_range_mask = (
            packets_df_line_interaction["nus"]
            < packet_nu_range[0]
        ) & (
            packets_df_line_interaction["nus"]
            > packet_nu_range[1]
        )

     # Calculate the area term to convert luminosity to flux
    if distance is None:
        lum_to_flux = 1  # so that this term will have no effect
    else:
        if distance <= 0:
            raise ValueError(
                "distance passed must be greater than 0. If you intended "
                "to plot luminosities instead of flux, set distance=None "
                "or don't specify distance parameter in the function call."
            )
        else:
            lum_to_flux = 4.0 * np.pi * (distance.to("cm")) ** 2


    # Histogram weights are packet luminosities or flux
    weights = (
        packets_df["energies"][packet_nu_range_mask]
        / lum_to_flux
    ) / time_of_simulation

    luminosities_df = pd.DataFrame(index=plot_wavelength)

    # Contribution of packets which experienced no interaction
    mask_noint = (
        packets_df["last_interaction_type"][packet_nu_range_mask]
        == -1
    )

    hist_noint = np.histogram(
        packets_df["nus"][packet_nu_range_mask][mask_noint],
        bins=plot_frequency_bins.value,
        weights=weights[mask_noint],
        density=False,
    )

    L_nu_noint = (
        hist_noint[0]
        * u.erg
        / u.s
        / spectrum_delta_frequency
    )
    L_lambda_noint = L_nu_noint * plot_frequency / plot_wavelength

    luminosities_df["noint"] = L_lambda_noint.value

    # Contribution of packets which only experienced electron scattering
    mask_escatter = (
        packets_df["last_interaction_type"][packet_nu_range_mask]
        == 1
    ) & (
        packets_df["last_line_interaction_in_id"][packet_nu_range_mask]
        == -1
    )
    hist_escatter = np.histogram(
        packets_df["nus"][packet_nu_range_mask][mask_escatter],
        bins=plot_frequency_bins.value,
        weights=weights[mask_escatter],
        density=False,
    )

    L_nu_escatter = (
        hist_escatter[0]
        * u.erg
        / u.s
        / spectrum_delta_frequency
    )
    L_lambda_escatter = L_nu_escatter * plot_frequency / plot_wavelength
    luminosities_df["escatter"] = L_lambda_escatter.value

    # Group packets_df by atomic number of elements with which packets
    # had their last emission (interaction out)
    if species_list is None:
        packets_df_grouped = (
            packets_df_line_interaction.loc[packet_nu_line_range_mask]
            .groupby(by="last_line_interaction_atom")
        )
    else:
        packets_df_grouped = (
            packets_df_line_interaction.loc[packet_nu_line_range_mask]
            .groupby(by="last_line_interaction_species")
        )

    # Contribution of each species with which packets interacted
    for identifier, group in packets_df_grouped:
        hist_el = np.histogram(
            group["nus"],
            bins=plot_frequency_bins.value,
            weights=group["energies"]
            / lum_to_flux
            / time_of_simulation,
        )

        L_nu_el = (
            hist_el[0]
            * u.erg
            / u.s
            / spectrum_delta_frequency
        )
        L_lambda_el = L_nu_el * plot_frequency / plot_wavelength

        luminosities_df[identifier] = L_lambda_el.value

    # Create an array of the species with which packets interacted
    emission_species_present = np.array(list(packets_df_grouped.groups.keys()))

    return luminosities_df, emission_species_present
