import logging

import astropy.units as u
import numpy as np
import pandas as pd

from tardis.energy_input.energy_source import (
    get_nuclear_lines_database,
)
from tardis.energy_input.gamma_packet_loop import gamma_packet_loop
from tardis.energy_input.gamma_ray_channel import (
    calculate_total_decays,
    create_inventories_dict,
    create_isotope_dicts,
)

from tardis.energy_input.gamma_ray_transport import (
    calculate_total_decays_old,
    create_isotope_dicts_old,
    create_inventories_dict_old,
)
from tardis.energy_input.gamma_ray_packet_source import RadioactivePacketSource
from tardis.energy_input.gamma_ray_transport import (
    calculate_average_energies,
    calculate_average_power_per_mass,
    calculate_ejecta_velocity_volume,
    calculate_energy_per_mass,
    decay_chain_energies,
    distribute_packets,
    get_taus,
    iron_group_fraction_per_shell,
)
from tardis.energy_input.GXPacket import GXPacket

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def run_gamma_ray_loop(
    model,
    plasma,
    num_decays,
    time_start,
    time_end,
    time_space,
    time_steps,
    seed,
    positronium_fraction,
    path_to_decay_data,
    spectrum_bins,
    grey_opacity,
    photoabsorption_opacity="tardis",
    pair_creation_opacity="tardis",
):
    """
    Main loop to determine the gamma-ray propagation through the ejecta.
    """
    np.random.seed(seed)
    time_explosion = model.time_explosion.to(u.s).value
    inner_velocities = model.v_inner.to("cm/s").value
    outer_velocities = model.v_outer.to("cm/s").value
    ejecta_volume = model.volume.to("cm^3").value
    number_of_shells = model.no_of_shells
    shell_masses = model.volume * model.density
    raw_isotope_abundance = model.composition.raw_isotope_abundance.sort_values(
        by=["atomic_number", "mass_number"], ascending=False
    )
    time_start *= u.d.to(u.s)
    time_end *= u.d.to(u.s)

    assert time_start < time_end, "time_start must be smaller than time_end!"
    if time_space == "log":
        times = np.geomspace(time_start, time_end, time_steps + 1)
        effective_time_array = np.sqrt(times[:-1] * times[1:])
    else:
        times = np.linspace(time_start, time_end, time_steps + 1)
        effective_time_array = 0.5 * (times[:-1] + times[1:])

    dt_array = np.diff(times)

    ejecta_velocity_volume = calculate_ejecta_velocity_volume(model)

    inv_volume_time = (
        1.0 / ejecta_velocity_volume[:, np.newaxis]
    ) / effective_time_array**3.0

    energy_df_rows = np.zeros((number_of_shells, time_steps))
    # Use isotopic number density
    for atom_number in plasma.isotope_number_density.index.get_level_values(0):
        values = plasma.isotope_number_density.loc[atom_number].values
        if values.shape[0] > 1:
            plasma.isotope_number_density.loc[atom_number].update = np.sum(
                values, axis=0
            )
        else:
            plasma.isotope_number_density.loc[atom_number].update = values

    # Electron number density
    electron_number_density = plasma.number_density.mul(
        plasma.number_density.index, axis=0
    ).sum()
    taus, parents = get_taus(raw_isotope_abundance)
    # inventories = raw_isotope_abundance.to_inventories()
    electron_number = np.array(electron_number_density * ejecta_volume)
    electron_number_density_time = (
        electron_number[:, np.newaxis] * inv_volume_time
    )

    # Calculate decay chain energies

    mass_density_time = shell_masses[:, np.newaxis] * inv_volume_time
    gamma_ray_lines = get_nuclear_lines_database(path_to_decay_data)
    isotope_dict = create_isotope_dicts_old(raw_isotope_abundance, shell_masses)
    inventories_dict = create_inventories_dict_old(isotope_dict)
    total_decays = calculate_total_decays_old(
        inventories_dict, time_end - time_start
    )

    (
        average_energies,
        average_positron_energies,
        gamma_ray_line_dict,
    ) = calculate_average_energies(raw_isotope_abundance, gamma_ray_lines)

    decayed_energy = decay_chain_energies(
        average_energies,
        total_decays,
    )
    energy_per_mass, energy_df = calculate_energy_per_mass(
        decayed_energy, raw_isotope_abundance, shell_masses
    )
    average_power_per_mass = calculate_average_power_per_mass(
        energy_per_mass, time_end - time_start
    )
    number_of_isotopes = plasma.isotope_number_density * ejecta_volume
    total_isotope_number = number_of_isotopes.sum().sum()
    decayed_packet_count = num_decays * number_of_isotopes.divide(
        total_isotope_number, axis=0
    )

    total_energy = energy_df.sum().sum()
    energy_per_packet = total_energy / num_decays
    packets_per_isotope_df = (
        distribute_packets(decayed_energy, total_energy, num_decays)
        .round()
        .fillna(0)
        .astype(int)
    )

    total_energy = total_energy * u.eV.to("erg")

    logger.info(f"Total gamma-ray energy is {total_energy}")

    iron_group_fraction = iron_group_fraction_per_shell(model)
    number_of_packets = packets_per_isotope_df.sum().sum()
    logger.info(f"Total number of packets is {number_of_packets}")
    individual_packet_energy = total_energy / number_of_packets
    logger.info(f"Energy per packet is {individual_packet_energy}")

    logger.info("Initializing packets")

    packet_source = RadioactivePacketSource(
        individual_packet_energy,
        gamma_ray_line_dict,
        positronium_fraction,
        inner_velocities,
        outer_velocities,
        inv_volume_time,
        times,
        energy_df_rows,
        effective_time_array,
        taus,
        parents,
        average_positron_energies,
        average_power_per_mass,
    )

    packet_collection = packet_source.create_packets(packets_per_isotope_df)

    energy_df_rows = packet_source.energy_df_rows
    energy_plot_df_rows = np.zeros((number_of_packets, 8))

    logger.info("Creating packet list")
    packets = []
    total_cmf_energy = packet_collection.energy_cmf.sum()
    total_rf_energy = packet_collection.energy_rf.sum()
    for i in range(number_of_packets):
        packet = GXPacket(
            packet_collection.location[:, i],
            packet_collection.direction[:, i],
            packet_collection.energy_rf[i],
            packet_collection.energy_cmf[i],
            packet_collection.nu_rf[i],
            packet_collection.nu_cmf[i],
            packet_collection.status[i],
            packet_collection.shell[i],
            packet_collection.time_current[i],
        )
        packets.append(packet)
        energy_plot_df_rows[i] = np.array(
            [
                i,
                packet.energy_rf,
                packet.get_location_r(),
                packet.time_current,
                int(packet.status),
                0,
                0,
                0,
            ]
        )

    logger.info(f"Total cmf energy is {total_cmf_energy}")
    logger.info(f"Total rf energy is {total_rf_energy}")

    energy_bins = np.logspace(2, 3.8, spectrum_bins)
    energy_out = np.zeros((len(energy_bins - 1), time_steps))
    packets_info_array = np.zeros((int(num_decays), 8))

    (
        energy_df_rows,
        energy_plot_df_rows,
        energy_out,
        deposition_estimator,
        bin_width,
        packets_array,
    ) = gamma_packet_loop(
        packets,
        grey_opacity,
        photoabsorption_opacity,
        pair_creation_opacity,
        electron_number_density_time,
        mass_density_time,
        inv_volume_time,
        iron_group_fraction.to_numpy(),
        inner_velocities,
        outer_velocities,
        times,
        dt_array,
        effective_time_array,
        energy_bins,
        energy_df_rows,
        energy_plot_df_rows,
        energy_out,
        packets_info_array,
    )

    energy_plot_df = pd.DataFrame(
        data=energy_plot_df_rows,
        columns=[
            "packet_index",
            "energy_input",
            "energy_input_r",
            "energy_input_time",
            "energy_input_type",
            "compton_opacity",
            "photoabsorption_opacity",
            "total_opacity",
        ],
    )

    energy_plot_positrons = pd.DataFrame(
        data=packet_source.energy_plot_positron_rows,
        columns=[
            "packet_index",
            "energy_input",
            "energy_input_r",
            "energy_input_time",
        ],
    )

    energy_estimated_deposition = (
        pd.DataFrame(data=deposition_estimator, columns=times[:-1])
    ) / dt_array

    energy_df = pd.DataFrame(data=energy_df_rows, columns=times[:-1]) / dt_array
    escape_energy = pd.DataFrame(
        data=energy_out, columns=times[:-1], index=energy_bins
    )

    return (
        energy_df,
        energy_plot_df,
        escape_energy,
        decayed_packet_count,
        energy_plot_positrons,
        energy_estimated_deposition,
    )
