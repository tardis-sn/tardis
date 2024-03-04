import numpy as np
import pandas as pd
import astropy.units as u

from tardis.energy_input.energy_source import (
    get_all_isotopes,
    get_nuclear_lines_database,
)

from tardis.energy_input.gamma_packet_loop import gamma_packet_loop
from tardis.energy_input.gamma_ray_transport import (
    initialize_packets,
    calculate_shell_masses,
    calculate_ejecta_velocity_volume,
    create_isotope_dicts,
    create_inventories_dict,
    calculate_total_decays,
    calculate_average_energies,
    decay_chain_energies,
    calculate_energy_per_mass,
    calculate_average_power_per_mass,
    iron_group_fraction_per_shell,
    get_taus,
    fractional_decay_energy,
    packets_per_isotope,
)


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
    number_of_shells = model.no_of_shells
    ejecta_volume = model.volume.to("cm^3").value
    raw_isotope_abundance = (
        model.composition.isotopic_mass_fraction.sort_values(
            by=["atomic_number", "mass_number"], ascending=False
        )
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
    inventories = raw_isotope_abundance.to_inventories()
    electron_number = np.array(electron_number_density * ejecta_volume)
    electron_number_density_time = (
        electron_number[:, np.newaxis] * inv_volume_time
    )

    # Calculate decay chain energies
    shell_masses = calculate_shell_masses(model)
    mass_density_time = shell_masses[:, np.newaxis] * inv_volume_time
    gamma_ray_lines = get_nuclear_lines_database(path_to_decay_data)
    isotope_dict = create_isotope_dicts(raw_isotope_abundance, shell_masses)
    inventories_dict = create_inventories_dict(isotope_dict)
    total_decays = calculate_total_decays(
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

    decayed_packet_count_dict = decayed_packet_count.to_dict()
    fractional_decay_energy_dict = fractional_decay_energy(decayed_energy)
    packets_per_isotope_df = (
        packets_per_isotope(
            fractional_decay_energy_dict, decayed_packet_count_dict
        )
        .round()
        .fillna(0)
        .astype(int)
    )
    print(packets_per_isotope_df)
    total_energy = energy_df.sum().sum()
    print("Total gamma-ray energy:")
    print(total_energy * u.keV.to("erg"))

    iron_group_fraction = iron_group_fraction_per_shell(model)
    number_of_packets = packets_per_isotope_df.sum().sum()
    print("Total number of packets:", number_of_packets)
    individual_packet_energy = total_energy / number_of_packets
    print("Energy per packet:", individual_packet_energy)

    print("Initializing packets")

    (
        packets,
        energy_df_rows,
        energy_plot_df_rows,
        energy_plot_positron_rows,
    ) = initialize_packets(
        packets_per_isotope_df,
        individual_packet_energy,
        positronium_fraction,
        inner_velocities,
        outer_velocities,
        gamma_ray_line_dict,
        average_positron_energies,
        inv_volume_time,
        times,
        energy_df_rows,
        effective_time_array,
        taus,
        parents,
        inventories,
        average_power_per_mass,
    )

    print("Total positron energy from packets")
    print(energy_df_rows.sum().sum() * u.eV.to("erg"))

    total_cmf_energy = 0.0
    total_rf_energy = 0.0

    for p in packets:
        total_cmf_energy += p.energy_cmf
        total_rf_energy += p.energy_rf

    print("Total cmf energy")
    print(total_cmf_energy)
    print("Total rf energy")
    print(total_rf_energy)

    energy_bins = np.logspace(2, 3.8, spectrum_bins)
    energy_out = np.zeros((len(energy_bins - 1), time_steps))
    packets_out = np.zeros((len(energy_bins - 1), 5))
    packets_info_array = np.zeros((num_decays, 10))

    (
        energy_df_rows,
        energy_plot_df_rows,
        energy_out,
        deposition_estimator,
        packets_out,
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
        packets_out,
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
        data=energy_plot_positron_rows,
        columns=[
            "packet_index",
            "energy_input",
            "energy_input_r",
            "energy_input_time",
        ],
    )

    packets_out_df = pd.DataFrame(
        data=packets_array,
        columns=[
            "number",
            "status",
            "Z",
            "A",
            "nu_cmf",
            "nu_rf",
            "energy_cmf",
            "lum_rf",
            "energy_rf",
            "shell",
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
        packets_out_df,
    )
