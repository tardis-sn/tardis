import logging
from typing import Optional

import astropy.units as u
import numpy as np
import pandas as pd
from numba import objmode, njit, prange

from numba.np.ufunc.parallel import get_num_threads, get_thread_id

from tardis.energy_input.gamma_ray_packet_source import (
    GammaRayPacketSource,
    legacy_calculate_positron_fraction,
)
from tardis.energy_input.gamma_ray_transport import (
    calculate_ejecta_velocity_volume,
    iron_group_fraction_per_shell,
)
from tardis.energy_input.transport.gamma_packet_loop import (
    single_gamma_packet_loop,
)
from tardis.energy_input.transport.GXPacket import GXPacket
from tardis.energy_input.util import get_index
from tardis.model.base import SimulationState
from tardis.transport.montecarlo import njit_dict, njit_dict_no_parallel

from tardis.energy_input.util import (
    C_CGS,
    H_CGS_KEV,
)

from tardis.util.base import update_packet_pbar

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def calculate_electron_number_density(
    simulation_state: SimulationState,
    ejecta_volume: np.ndarray,
    effective_time_array: np.ndarray,
    legacy: bool = False,
    legacy_atom_data: Optional[object] = None,
) -> np.ndarray:
    """
    Calculate electron number density and its time evolution.

    Parameters
    ----------
    model : SimulationState
        Tardis model object.
    ejecta_volume : np.ndarray
        Ejecta volume array.
    effective_time_array : np.ndarray
        Effective time array in seconds.
    legacy : bool, optional
        Whether to use legacy calculation method. Default is False.
    legacy_atom_data : object, optional
        Legacy atom data object. Required if legacy=True.

    Returns
    -------
    electron_number_density_time : np.ndarray
        Electron number density evolution with time.
    """
    ejecta_velocity_volume = calculate_ejecta_velocity_volume(simulation_state)

    inv_volume_time = (
        1.0 / ejecta_velocity_volume[:, np.newaxis]
    ) / effective_time_array**3.0

    # Calculate the elemental number density
    if not legacy:
        elemental_number_density = (
            simulation_state.composition.isotopic_number_density.groupby(
                "atomic_number"
            ).sum()
        )
    else:
        if legacy_atom_data is None:
            raise ValueError("legacy_atom_data must be provided when legacy=True")
        elemental_number_density = simulation_state.calculate_elemental_number_density(
            legacy_atom_data.atom_data.mass
        )

    # Electron number density
    electron_number_density = elemental_number_density.mul(
        elemental_number_density.index,
        axis=0,
    ).sum()
    electron_number = np.array(electron_number_density * ejecta_volume)

    # Evolve electron number and mass density with time
    electron_number_density_time = electron_number[:, np.newaxis] * inv_volume_time

    return electron_number_density_time


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
    simulation_state,
    legacy_isotope_decacy_df,
    cumulative_decays_df,
    number_of_packets,
    times,
    effective_time_array,
    seed,
    positronium_fraction,
    spectrum_bins,
    grey_opacity,
    photo_absorption_opacity="tardis",
    pair_creation_opacity="tardis",
    legacy=False,
    legacy_atom_data=None,
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
    inner_velocities = simulation_state.v_inner.to("cm/s").value
    outer_velocities = simulation_state.v_outer.to("cm/s").value
    ejecta_volume = simulation_state.volume.to("cm^3").value
    shell_masses = simulation_state.volume * simulation_state.density
    number_of_shells = len(shell_masses)
    # TODO: decaying upto times[0]. raw_isotope_abundance is possibly not the best name
    isotopic_mass_fraction = (
        simulation_state.composition.isotopic_mass_fraction.sort_values(
            by=["atomic_number", "mass_number"], ascending=False
        )
    )

    dt_array = np.diff(times)

    # Calculate electron number density evolution
    electron_number_density_time = calculate_electron_number_density(
        simulation_state,
        ejecta_volume,
        effective_time_array,
        legacy=legacy,
        legacy_atom_data=legacy_atom_data,
    )

    # Calculate mass density evolution
    ejecta_velocity_volume = calculate_ejecta_velocity_volume(simulation_state)
    inv_volume_time = (
        1.0 / ejecta_velocity_volume[:, np.newaxis]
    ) / effective_time_array**3.0
    mass_density_time = shell_masses[:, np.newaxis] * inv_volume_time

    packet_source = GammaRayPacketSource(
        cumulative_decays_df=cumulative_decays_df,
        isotope_decay_df=legacy_isotope_decacy_df,
        positronium_fraction=positronium_fraction,
        inner_velocities=inner_velocities,
        outer_velocities=outer_velocities,
        times=times,
        effective_times=effective_time_array,
        base_seed=seed,
    )

    logger.info("Creating packets")
    if legacy:
        # Calculate energy per packet for legacy mode using legacy_isotope_decacy_df
        gamma_df = legacy_isotope_decacy_df[
            legacy_isotope_decacy_df["radiation"] == "g"
        ]
        total_energy_gamma = gamma_df["decay_energy_erg"].sum()
        energy_per_packet = total_energy_gamma / number_of_packets
        legacy_energy_per_packet = energy_per_packet
    else:
        # Let the packet source calculate energy per packet internally
        legacy_energy_per_packet = None
        energy_per_packet = None

    packet_collection, source_isotopes = packet_source.create_packets(
        cumulative_decays_df,
        number_of_packets,
        legacy_energy_per_packet=legacy_energy_per_packet,
    )

    # For non-legacy mode, get the energy per packet from the packet source calculation
    if not legacy:
        gamma_df = cumulative_decays_df[cumulative_decays_df["radiation"] == "g"]
        total_energy_gamma = gamma_df["decay_energy_erg"].sum()
        energy_per_packet = total_energy_gamma / number_of_packets

    total_energy = np.zeros((number_of_shells, len(times) - 1))

    #logger.info("Creating packet list")
    #packets = []
    ## This for loop is expensive. Need to rewrite GX packet to handle arrays

#    #packets = get_packet_properties_parallel(
    #    packet_collection, number_of_packets)

#    ## Calculate isotope positron fraction separately
    #isotope_positron_fraction = legacy_calculate_positron_fraction(
    #    legacy_isotope_decacy_df, packet_collection.source_isotopes, number_of_packets
    #)
    #for i, p in enumerate(packets):
    #    total_energy[p.shell, p.time_index] += isotope_positron_fraction[i] * energy_per_packet

#    #logger.info(
    #    "Total energy deposited by the positrons is %s", total_energy.sum().sum()
    # )

    # positron_energy = total_energy

    # Copy the positron energy to a new dataframe
    #positron_energy_df = pd.DataFrame(data=total_energy.copy(), columns=times[:-1])

    number_of_threads = get_num_threads()

    #for packet_index in range(number_of_packets):
    energy_bins = np.logspace(2, 3.8, spectrum_bins)
    energy_out = np.zeros((number_of_threads, len(energy_bins) - 1, len(times) - 1))
    energy_out_cosi = np.zeros((len(energy_bins) - 1, len(times) - 1))
    energy_deposited = np.zeros((number_of_threads, number_of_shells , len(times) - 1))
    packets_info_array = np.zeros((int(number_of_packets), 8))
    iron_group_fraction = iron_group_fraction_per_shell(simulation_state).to_numpy()

    logger.info("Entering the main gamma-ray loop")

    total_cmf_energy = 0
    total_rf_energy = 0

    (energy_out,
    energy_deposited_gamma) = run_gamma_ray_loop_parallel(packet_collection,
                                number_of_packets,
                                grey_opacity,
                                photo_absorption_opacity,
                                pair_creation_opacity,
                                electron_number_density_time,
                                mass_density_time,
                                iron_group_fraction,
                                inner_velocities,
                                outer_velocities,
                                dt_array,
                                times,
                                effective_time_array,
                                energy_bins,
                                energy_out,
                                energy_out_cosi,
                                total_energy,
                                energy_deposited,
                                )
    return (energy_out,
            energy_deposited_gamma
           )



    #for p in packets:
    #    total_cmf_energy += p.energy_cmf
    #    total_rf_energy += p.energy_rf

    #logger.info("Total CMF energy is %s", total_cmf_energy)
    #logger.info("Total RF energy is %s", total_rf_energy)

    #(
    #    energy_out,
    #    energy_out_cosi

    # Run the gamma-ray loop in parallel    (

    #    packets_array,
    #    energy_deposited_gamma,
    #    total_energy,
    #) = gamma_packet_loop(
    #    packets,
    #    grey_opacity,
    #    photoabsorption_opacity,
    #    pair_creation_opacity,
    #     electron_number_density_time,
    #     mass_density_time,
    #     iron_group_fraction.to_numpy(),
    #     inner_velocities,
    #     outer_velocities,
    #     dt_array,
    #     times,
    #     effective_time_array,
    #     energy_bins,
    #     energy_out,
    #     energy_out_cosi,
    #     total_energy,
    #     energy_deposited,
    #     packets_info_array,
    # )

#     # packets_df_escaped = pd.DataFrame(
    #     data=packets_array,
    #     columns=[
    #         "packet_index",
    #         "status",
    #         "nu_cmf",
    #         "nu_rf",
    #         "energy_cmf",
    #         "luminosity",
    #         "energy_rf",
    #        "shell_number",
    #     ],
    # )

#     # escape_energy = pd.DataFrame(
    #     data=energy_out, columns=times[:-1], index=energy_bins
    # )

#     # escape_energy_cosi = pd.DataFrame(
    #     data=energy_out_cosi, columns=times[:-1], index=energy_bins
    # )

#     # # deposited energy by gamma-rays in ergs
    # gamma_ray_deposited_energy = pd.DataFrame(
    #    data=energy_deposited_gamma, columns=times[:-1]
    # )
    # # deposited energy by positrons and gamma-rays in ergs
    # total_energy = pd.DataFrame(
    #    data=total_energy, columns=times[:-1]
    # )

#     # logger.info(
    #    "Total energy deposited by gamma rays and positrons is %s",
    #     total_energy.sum().sum(),
    # )

#     # total_deposited_energy = total_energy / dt_array

#     # return (
    #     escape_energy,
    #     escape_energy_cosi,
    #    packets_df_escaped,
    #     gamma_ray_deposited_energy,
    #     total_deposited_energy,
    #     positron_energy_df
    # )


# def get_packet_properties(number_of_shells, times, time_steps, packets):
    # """
    # Function to get the properties of the packets.

#     # Parameters
    # ----------
    # packets : list
    #     List of packets.


# 

@njit(**njit_dict_no_parallel)
def run_gamma_ray_loop_parallel(packet_collection,
                                number_of_packets,
                                grey_opacity,
                                photo_absorption_opacity,
                                pair_creation_opacity,
                                electron_number_density_time,
                                mass_density_time,
                                iron_group_fraction,
                                inner_velocities,
                                outer_velocities,
                                dt_array,
                                times,
                                effective_time_array,
                                energy_bins,
                                energy_out,
                                energy_out_cosi,
                                total_energy,
                                energy_deposited_gamma,
                                show_progress_bars=True
                                ):
    
    """    
    Function to run the gamma-ray loop in parallel.
    """

    main_thread_id = get_thread_id()
    n_threads = get_num_threads()

    for i in prange(number_of_packets):
        thread_id = get_thread_id()
        if show_progress_bars:
            if thread_id == main_thread_id:
                with objmode:
                    update_amount = 1 * n_threads
                    update_packet_pbar(
                        update_amount,
                        current_iteration=1,
                        no_of_packets=number_of_packets,
                        total_iterations=1,
                    )


        gx_packet = GXPacket(
           packet_collection.location[:, i],
           packet_collection.direction[:, i],
           packet_collection.energy_rf[i],
           packet_collection.energy_cmf[i],
           packet_collection.nu_rf[i],
           packet_collection.nu_cmf[i],
           packet_collection.status[i],
           packet_collection.shell[i],
           packet_collection.time_start[i],
           packet_collection.time_index[i],
       )
        
        (
        packet_info, 
        ejecta_energy_gained, 
        ) = single_gamma_packet_loop(
            gx_packet,
            grey_opacity,
            photo_absorption_opacity,
            pair_creation_opacity,
            electron_number_density_time,
            mass_density_time,
            iron_group_fraction,
            inner_velocities,
            outer_velocities,
            dt_array,
            times,
            effective_time_array,
            energy_bins,
            energy_out,
            energy_out_cosi,
            total_energy,
            energy_deposited_gamma,
        )
        # print (gx_packet.shell)
        if gx_packet.status != 5:
            energy_deposited_gamma[thread_id, gx_packet.shell, gx_packet.time_index] += ejecta_energy_gained
        else:
            pass



        rest_energy = gx_packet.nu_rf * H_CGS_KEV
        bin_index = get_index(rest_energy, energy_bins)
        bin_width = energy_bins[bin_index + 1] - energy_bins[bin_index]
        freq_bin_width = bin_width / H_CGS_KEV
        energy_out[thread_id, bin_index, gx_packet.time_index] += (
            gx_packet.energy_rf
                    / dt_array[gx_packet.time_index]
                    / freq_bin_width
                )
        #total_energy[i, gx_packet.shell, time_index] = ejecta_energy_gained
    
    energy_deposited_gamma = energy_deposited_gamma.sum(axis=0)
    energy_out = energy_out.sum(axis=0)
    # Fix with positrons energy later
    #total_energy = energy_deposited_gamma

    return (energy_out,
        energy_deposited_gamma,
        )