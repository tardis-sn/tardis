import logging

import astropy.units as u
import numpy as np
import pandas as pd

from tardis.energy_input.energy_source import (
    get_nuclear_lines_database,
)
from tardis.energy_input.gamma_packet_loop import gamma_packet_loop


from tardis.energy_input.gamma_ray_transport import (
    calculate_total_decays_old,
    create_isotope_dicts_old,
    create_inventories_dict_old,
)
from tardis.energy_input.gamma_ray_packet_source import (
    RadioactivePacketSource,
    GammaRayPacketSource,
)
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
from tardis.energy_input.util import make_isotope_string_tardis_like

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def get_effective_time_array(time_start, time_end, time_space, time_steps):
    """
    Function to get the effective time array for the gamma-ray loop.

    Parameters
    ----------
    time_start : float
        start time in days.
    time_end : float
        end time in days.
    time_space : str
        linear or log.
    time_steps : int
        number of time steps.

    Returns
    -------
    times : np.ndarray
        array of times in secs.
    effective_time_array : np.ndarray
        effective time array in secs.
    """

    assert time_start < time_end, "time_start must be smaller than time_end!"
    if time_space == "log":
        times = np.geomspace(time_start, time_end, time_steps + 1)
        effective_time_array = np.sqrt(times[:-1] * times[1:])
    else:
        times = np.linspace(time_start, time_end, time_steps + 1)
        effective_time_array = 0.5 * (times[:-1] + times[1:])

    return times, effective_time_array


def run_gamma_ray_loop(
    model,
    plasma,
    isotope_decay_df,
    cumulative_decays_df,
    num_decays,
    times,
    effective_time_array,
    seed,
    positronium_fraction,
    spectrum_bins,
    time_steps,
    grey_opacity,
    photoabsorption_opacity="tardis",
    pair_creation_opacity="tardis",
):
    """
    Main loop to determine the gamma-ray propagation through the ejecta.

    Parameters
    ----------
    model : tardis.model.Radial1DModel
            Tardis model object.
    plasma : tardis.plasma.standard_plasmas.BasePlasma
            Tardis plasma object.
    isotope_decay_df : pd.DataFrame
            DataFrame containing the cumulative decay data.
    cumulative_decays_df : pd.DataFrame
            DataFrame containing the time evolving mass fractions.
    num_decays : int
            Number of packets to decay.
    times : np.ndarray
            Array of times in days.
    effective_time_array : np.ndarray
            Effective time array in days.
    seed : int
            Seed for the random number generator.
    positronium_fraction : float
            Fraction of positronium.
    spectrum_bins : int
            Number of spectrum bins.
    time_steps : int
            Number of time steps.
    grey_opacity : float
            Grey opacity.
    photoabsorption_opacity : str
            Photoabsorption opacity.
    pair_creation_opacity : str
            Pair creation opacity.


    Returns
    -------

    escape_energy : pd.DataFrame
            DataFrame containing the energy escaping the ejecta.
    packets_df_escaped : pd.DataFrame
            DataFrame containing the packets info that escaped the ejecta.


    """
    np.random.seed(seed)
    times = times * u.d.to(u.s)
    effective_time_array = effective_time_array * u.d.to(u.s)
    inner_velocities = model.v_inner.to("cm/s").value
    outer_velocities = model.v_outer.to("cm/s").value
    ejecta_volume = model.volume.to("cm^3").value
    shell_masses = model.volume * model.density
    raw_isotope_abundance = model.composition.raw_isotope_abundance.sort_values(
        by=["atomic_number", "mass_number"], ascending=False
    )

    dt_array = np.diff(times)

    ejecta_velocity_volume = calculate_ejecta_velocity_volume(model)

    inv_volume_time = (
        1.0 / ejecta_velocity_volume[:, np.newaxis]
    ) / effective_time_array**3.0

    for (
        atom_number
    ) in model.composition.isotopic_number_density.index.get_level_values(0):
        values = model.composition.isotopic_number_density.loc[
            atom_number
        ].values
        if values.shape[0] > 1:
            plasma.number_density.loc[atom_number].update = np.sum(
                values, axis=0
            )
        else:
            plasma.number_density.loc[atom_number].update = values

    # Electron number density
    electron_number_density = plasma.number_density.mul(
        plasma.number_density.index,
        axis=0,
    ).sum()
    electron_number = np.array(electron_number_density * ejecta_volume)

    # Evolve electron number and mass density with time
    electron_number_density_time = (
        electron_number[:, np.newaxis] * inv_volume_time
    )
    mass_density_time = shell_masses[:, np.newaxis] * inv_volume_time

    taus, parents = get_taus(raw_isotope_abundance)
    # Need to get the strings for the isotopes without the dashes
    taus = make_isotope_string_tardis_like(taus)

    # Get only the gamma rays
    isotope_decay_df_gamma = isotope_decay_df[
        isotope_decay_df["radiation"] == "g"
    ]
    total_energy_gamma = isotope_decay_df_gamma["decay_energy_erg"].sum()

    energy_per_packet = total_energy_gamma / num_decays

    logger.info(f"Total energy in gamma-rays is {total_energy_gamma}")
    logger.info(f"Energy per packet is {energy_per_packet}")

    packet_source = GammaRayPacketSource(
        energy_per_packet,
        isotope_decay_df,
        positronium_fraction,
        inner_velocities,
        outer_velocities,
        inv_volume_time,
        times,
        effective_time_array,
        taus,
        parents,
    )
    logger.info("Creating packets")
    packet_collection = packet_source.create_packets(
        cumulative_decays_df, num_decays, seed
    )
    logger.info("Creating packet list")
    packets = []
    packets = [
        GXPacket(
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
        for i in range(num_decays)
    ]

    energy_bins = np.logspace(2, 3.8, spectrum_bins)
    energy_out = np.zeros((len(energy_bins - 1), time_steps))
    packets_info_array = np.zeros((int(num_decays), 8))
    iron_group_fraction = iron_group_fraction_per_shell(model)

    logger.info("Entering the main gamma-ray loop")

    total_cmf_energy = 0
    total_rf_energy = 0

    for p in packets:
        total_cmf_energy += p.energy_cmf
        total_rf_energy += p.energy_rf

    logger.info(f"Total CMF energy is {total_cmf_energy}")
    logger.info(f"Total RF energy is {total_rf_energy}")

    (
        energy_out,
        packets_array,
    ) = gamma_packet_loop(
        packets,
        grey_opacity,
        photoabsorption_opacity,
        pair_creation_opacity,
        electron_number_density_time,
        mass_density_time,
        iron_group_fraction.to_numpy(),
        inner_velocities,
        outer_velocities,
        times,
        dt_array,
        effective_time_array,
        energy_bins,
        energy_out,
        packets_info_array,
    )

    packets_df_escaped = pd.DataFrame(
        data=packets_array,
        columns=[
            "packet_index",
            "status",
            "nu_cmf",
            "nu_rf",
            "energy_cmf",
            "luminosity",
            "energy_rf",
            "shell_number",
        ],
    )

    escape_energy = pd.DataFrame(
        data=energy_out, columns=effective_time_array, index=energy_bins
    )

    return escape_energy, packets_df_escaped
