import numpy as np
import pandas as pd
import astropy.units as u
from numba import njit
from numba.typed import List
import radioactivedecay as rd

from tardis.montecarlo.montecarlo_numba import njit_dict_no_parallel
from tardis.montecarlo.montecarlo_numba.opacities import M_P

from tardis.energy_input.energy_source import (
    get_all_isotopes,
    get_nuclear_lines_database,
    positronium_continuum,
    read_artis_lines,
    setup_input_energy,
)
from tardis.energy_input.samplers import initial_packet_radius
from tardis.energy_input.GXPacket import initialize_packet_properties
from tardis.energy_input.gamma_packet_loop import gamma_packet_loop

# Energy: keV, exported as eV for SF solver
# distance: cm
# mass: g
# time: s


def get_nuclide_atomic_number(nuclide):
    return rd.Nuclide(nuclide).Z


def get_chain_decay_power_per_ejectamass(
    inventory,
    time,
    time_start,
    average_energies,
    average_positron_energies,
    tau_end,
):
    # total decay power per mass [erg/s/g] for a given decaypath
    # only decays at the end of the chain contributed from the initial abundance of the top of the chain are counted
    # (these can be can be same for a chain of length one)

    # need to find a way to connect inventory_start with inventory_end to obtain start and end of chains
    # hopefully can just do it all at once

    # start_isotope is e.g. Ni56 while end_isotope would be e.g. Fe56
    inventory_end = inventory.decay(time - time_start)
    # contribution to the end nuclide abundance from the top of chain (could be a length-one chain Z,A_top = Z,A_end
    # so contribution would be from init abundance only)

    decaypower = []
    progeny = inventory.progeny()

    for start_isotope in progeny:
        # ignore branching probs for now
        end_isotope = progeny[start_isotope][0]

        if rd.Nuclide(end_isotope).half_life("readable") == "stable":
            end_isotope = start_isotope

        # i.e. need to get end_isotope(s) by finding out which ones have no progeny left
        # this should be the total number of "end chain" isotopes after (time - time_start)
        endnucabund = inventory_end.mass_fractions()[end_isotope]

        # this is the total decay energy from gamma-rays and positrons for the end of chain isotope
        # endecay = get_decaypath_lastnucdecayenergy(decaypathindex)
        endecay = (
            average_energies[end_isotope]
            + average_positron_energies[end_isotope]
        )

        print("Decay energy, abundance, tau")
        print(endecay)
        print(endnucabund)
        print(tau_end)

        decaypower.append(
            endecay
            * endnucabund
            / tau_end
            / (rd.Nuclide(start_isotope).atomic_mass * M_P)
        )

    return decaypower


@njit(**njit_dict_no_parallel)
def calculate_positron_fraction(
    positron_energy, isotope_energy, isotope_intensity
):
    """Calculate the fraction of energy that an isotope
    releases as positron kinetic energy

    Parameters
    ----------
    positron_energy : float
        Average kinetic energy of positrons from decay
    isotope_energy : numpy array
        Photon energies released by the isotope
    isotope_intensity : numpy array
        Intensity of photon energy release

    Returns
    -------
    float
        Fraction of energy released as positron kinetic energy
    """
    return positron_energy / np.sum(isotope_energy * isotope_intensity)


def initialize_packets(
    decays_per_isotope,
    packet_energy,
    gamma_ray_lines,
    positronium_fraction,
    inner_velocities,
    outer_velocities,
    inv_volume_time,
    times,
    energy_df_rows,
    effective_times,
    taus,
    parents,
    average_positron_energies,
    inventories,
    average_power_per_mass,
):
    """Initialize a list of GXPacket objects for the simulation
    to operate on.

    Parameters
    ----------
    decays_per_isotope : array int64
        Number of decays per simulation shell per isotope
    input_energy : float64
        Total input energy from decay
    ni56_lines : array float64
        Lines and intensities for Ni56
    co56_lines : array float64
        Lines and intensities for Co56
    inner_velocities : array float64
        Inner velocities of the shells
    outer_velocities : array float64
        Outer velocities of the shells
    inv_volume_time : array float64
        Inverse volume with time
    times : array float64
        Simulation time steps
    energy_df_rows : list
        Setup list for energy DataFrame output
    effective_times : array float64
        Middle time of the time step
    taus : array float64
        Mean lifetime for each isotope

    Returns
    -------
    list
        List of GXPacket objects
    array
        Array of main output dataframe rows
    array
        Array of plotting output dataframe rows
    array
        Array of positron output dataframe rows
    """

    packets = List()

    number_of_packets = decays_per_isotope.sum().sum()
    decays_per_shell = decays_per_isotope.T.sum().values

    energy_plot_df_rows = np.zeros((number_of_packets, 8))
    energy_plot_positron_rows = np.zeros((number_of_packets, 4))

    positronium_energy, positronium_intensity = positronium_continuum()

    packet_index = 0
    for k, shell in enumerate(decays_per_shell):

        initial_radii = initial_packet_radius(
            shell, inner_velocities[k], outer_velocities[k]
        )

        isotope_packet_count_df = decays_per_isotope.iloc[k]

        i = 0
        for (
            isotope_name,
            isotope_packet_count,
        ) in isotope_packet_count_df.items():
            isotope_energy = gamma_ray_lines[isotope_name][0, :]
            isotope_intensity = gamma_ray_lines[isotope_name][1, :]
            isotope_positron_fraction = calculate_positron_fraction(
                average_positron_energies[isotope_name],
                isotope_energy,
                isotope_intensity,
            )
            tau_start = taus[isotope_name]

            if isotope_name in parents:
                tau_end = taus[parents[isotope_name]]
            else:
                tau_end = 0

            for c in range(isotope_packet_count):
                packet, decay_time_index = initialize_packet_properties(
                    isotope_energy,
                    isotope_intensity,
                    positronium_energy,
                    positronium_intensity,
                    positronium_fraction,
                    packet_energy,
                    k,
                    tau_start,
                    tau_end,
                    initial_radii[i],
                    times,
                    effective_times,
                    inventories[k],
                    average_power_per_mass,
                )

                energy_df_rows[k, decay_time_index] += (
                    isotope_positron_fraction * packet_energy * 1000
                )

                energy_plot_df_rows[packet_index] = np.array(
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

                energy_plot_positron_rows[packet_index] = [
                    packet_index,
                    isotope_positron_fraction * packet_energy * 1000,
                    # * inv_volume_time[packet.shell, decay_time_index],
                    packet.get_location_r(),
                    packet.time_current,
                ]

                packets.append(packet)

                i += 1
                packet_index += 1

    return (
        packets,
        energy_df_rows,
        energy_plot_df_rows,
        energy_plot_positron_rows,
    )


def main_gamma_ray_loop(
    num_decays,
    model,
    plasma,
    time_steps=10,
    time_end=80.0,
    grey_opacity=-1,
    spectrum_bins=500,
    time_space="log",
    photoabsorption_opacity="tardis",
    pair_creation_opacity="tardis",
    seed=1,
    path_to_decay_data="~/Downloads/tardisnuclear/decay_radiation.h5",
    positronium_fraction=0.0,
):
    """Main loop that determines the gamma ray propagation

    Parameters
    ----------
    num_decays : int
        Number of decays requested
    model : tardis.Radial1DModel
        The tardis model to calculate gamma ray propagation through
    plasma : tardis.plasma.BasePlasma
        The tardis plasma with calculated atomic number density
    time_steps : int
        Number of time steps requested
    time_end : float
        End time of simulation in days
    grey_opacity : float
        Grey photoabsorption opacity for gamma-rays in cm^2 g^-1, set to -1 to turn off
    spectrum_bins : int
        Number of energy bins for the gamma-ray spectrum
    time_space : str
        `'log'` for log-space time steps, otherwise time steps will be linear
    photoabsorption_opacity : str
        Set the photoabsorption opacity calculation.
        Defaults to Ambwani & Sutherland (1988) approximation.
        `'kasen'` uses the Kasen et al. 2006 method.
    pair_creation_opacity : str
        Set the pair creation opacity calculation. Defaults to Ambwani & Sutherland (1988) approximation.
        `'artis'` uses the ARTIS implementation of the Ambwani & Sutherland (1988) approximation.
    seed : int
        Sets the seed for the random number generator. Uses deprecated methods.
    path_to_decay_data : str
        The path to a decay radiation file from the `nuclear` package.

    Returns
    -------
    pandas.DataFrame
        Energy per shell per time in units of eV/s/cm^-3
    pandas.DataFrame
        Columns:
        packet index,
        Energy input per packet,
        radius of deposition,
        time of deposition,
        compton opacity,
        photoabsorption opacity,
        pair creation opacity
    pandas.DataFrame
        Energy of escaping packets
    numpy.ndarray
        Packets emitted per shell
    pandas.DataFrame
        Energy from positrons
    pandas.DataFrame
        Estimated energy deposition in units of keV/s/cm^-3
    """
    # Note: not best numpy practice, but works better in numba than the alternatives
    np.random.seed(seed)

    # Enforce cgs
    outer_velocities = model.v_outer.to("cm/s").value
    inner_velocities = model.v_inner.to("cm/s").value
    ejecta_density = model.density.to("g/cm^3").value
    ejecta_volume = model.volume.to("cm^3").value
    ejecta_velocity_volume = (
        4 * np.pi / 3 * (outer_velocities**3.0 - inner_velocities**3.0)
    )
    time_explosion = model.time_explosion.to("s").value
    number_of_shells = model.no_of_shells
    raw_isotope_abundance = model.raw_isotope_abundance.sort_values(
        by=["atomic_number", "mass_number"], ascending=False
    )

    shell_masses = ejecta_volume * ejecta_density

    time_start = time_explosion
    time_end *= u.d.to(u.s)

    assert (
        time_start < time_end
    ), "Error, simulation start time greater than end time!"

    if time_space == "log":
        times = np.zeros(time_steps + 1)

        # log time steps
        for i in range(time_steps + 1):
            times[i] = (
                np.log(time_start)
                + (np.log(time_end) - np.log(time_start)) / time_steps * i
            )
            times[i] = np.exp(times[i])
    else:
        times = np.linspace(time_start, time_end, time_steps + 1)

    dt_array = np.diff(times)
    effective_time_array = np.array(
        [np.sqrt(times[i] * times[i + 1]) for i in range(time_steps)]
    )

    # Use isotopic number density
    for atom_number in plasma.isotope_number_density.index.get_level_values(0):
        values = plasma.isotope_number_density.loc[atom_number].values
        if values.shape[1] > 1:
            plasma.number_density.loc[atom_number] = np.sum(values, axis=0)
        else:
            plasma.number_density.loc[atom_number] = values

    # Calculate electron number density
    electron_number_density = (
        plasma.number_density.mul(plasma.number_density.index, axis=0)
    ).sum()

    electron_number_density_time = np.zeros(
        (len(ejecta_velocity_volume), len(effective_time_array))
    )

    mass_density_time = np.zeros(
        (len(ejecta_velocity_volume), len(effective_time_array))
    )

    electron_number = (electron_number_density * ejecta_volume).to_numpy()

    inv_volume_time = np.zeros(
        (len(ejecta_velocity_volume), len(effective_time_array))
    )

    # Pre-calculate quantities as they change with time
    for i, t in enumerate(effective_time_array):
        inv_volume_time[:, i] = (1.0 / ejecta_velocity_volume) / (t**3.0)
        mass_density_time[:, i] = shell_masses * inv_volume_time[:, i]
        electron_number_density_time[:, i] = (
            electron_number * inv_volume_time[:, i]
        )

    energy_df_rows = np.zeros((number_of_shells, time_steps))

    # Calculate number of packets per shell based on the mass of isotopes
    number_of_isotopes = plasma.isotope_number_density * ejecta_volume
    total_number_isotopes = number_of_isotopes.sum(axis=1)

    inventories = raw_isotope_abundance.to_inventories()
    all_isotope_names = get_all_isotopes(raw_isotope_abundance)
    all_isotope_names.sort()

    gamma_ray_lines = get_nuclear_lines_database(path_to_decay_data)

    taus = {}
    parents = {}
    gamma_ray_line_array_list = []
    average_energies_list = []
    average_positron_energies_list = []

    for i, isotope in enumerate(all_isotope_names):
        nuclide = rd.Nuclide(isotope)
        taus[isotope] = nuclide.half_life() / np.log(2)
        child = nuclide.progeny()
        if child is not None:
            for c in child:
                if rd.Nuclide(c).half_life("readable") != "stable":
                    parents[c] = isotope

        energy, intensity = setup_input_energy(
            gamma_ray_lines[
                gamma_ray_lines.Isotope == isotope.replace("-", "")
            ],
            "g",
        )
        gamma_ray_line_array_list.append(np.stack([energy, intensity]))
        average_energies_list.append(np.sum(energy * intensity))
        positron_energy, positron_intensity = setup_input_energy(
            gamma_ray_lines[
                gamma_ray_lines.Isotope == isotope.replace("-", "")
            ],
            "bp",
        )
        average_positron_energies_list.append(
            np.sum(positron_energy * positron_intensity)
        )

    # Construct Numba typed dicts
    gamma_ray_line_arrays = {}
    average_energies = {}
    average_positron_energies = {}

    for iso, lines in zip(all_isotope_names, gamma_ray_line_array_list):
        gamma_ray_line_arrays[iso] = lines

    for iso, energy, positron_energy in zip(
        all_isotope_names, average_energies_list, average_positron_energies_list
    ):
        average_energies[iso] = energy
        average_positron_energies[iso] = positron_energy

    # urilight chooses to have 0 as the baseline for this calculation
    # but time_start may also be valid in which case decay time is time_end - time_start
    total_energy_list = []

    for shell, inv in enumerate(inventories):
        decayed_energy = {}
        total_decays = inv.cumulative_decays(time_end)
        for nuclide in total_decays:
            if nuclide in parents and nuclide != "Co-56" and nuclide != "Co-57":
                parent = parents[nuclide]
                if parent in parents:
                    parent = parents[parent]
                decayed_energy[parent] += (
                    total_decays[nuclide]
                    * average_energies[nuclide]
                    * shell_masses[shell]
                )
            else:
                decayed_energy[nuclide] = (
                    total_decays[nuclide]
                    * average_energies[nuclide]
                    * shell_masses[shell]
                )

        total_energy_list.append(decayed_energy)

    total_energy = pd.DataFrame(total_energy_list)

    total_energy_columns = total_energy.columns.to_list()

    total_energy = total_energy[
        sorted(
            total_energy_columns, key=get_nuclide_atomic_number, reverse=True
        )
    ]

    energy_per_mass = total_energy.divide(
        (raw_isotope_abundance * shell_masses).T.to_numpy(),
        axis=0,
    )

    # Time averaged energy per mass for constant packet count
    average_power_per_mass = energy_per_mass / (time_end - time_start)

    energy_per_mass_norm = energy_per_mass.divide(
        energy_per_mass.sum(axis=1), axis=0
    )  # .cumsum(axis=1)

    decayed_packet_count = num_decays * number_of_isotopes.divide(
        total_number_isotopes, axis=0
    )

    packets_per_isotope = (
        (energy_per_mass_norm * decayed_packet_count.T.values)
        .round()
        .fillna(0)
        .astype(int)
    )

    print("Total gamma-ray energy")
    print(total_energy.sum().sum() * u.keV.to("erg"))

    print("Total positron energy")
    print(total_energy["Co-56"].sum(axis=0) * 0.0337 * u.keV.to("erg"))

    # Taking iron group to be elements 21-30
    # Used as part of the approximations for photoabsorption and pair creation
    # Dependent on atomic data
    iron_group_fraction_per_shell = model.abundance.loc[(21):(30)].sum(axis=0)

    number_of_packets = packets_per_isotope.sum().sum()
    print("Total packets:", number_of_packets)

    packet_energy = total_energy.sum().sum() / number_of_packets

    print("Energy per packet", packet_energy)

    # Need to update volume for positron deposition to be time-dependent
    print("Initializing packets")
    (
        packets,
        energy_df_rows,
        energy_plot_df_rows,
        energy_plot_positron_rows,
    ) = initialize_packets(
        packets_per_isotope,
        packet_energy,
        gamma_ray_line_arrays,
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
        inventories,
        average_power_per_mass,
    )

    print("Total positron energy from packets")
    print((energy_df_rows).sum().sum() * u.eV.to("erg"))

    total_cmf_energy = 0
    total_rf_energy = 0

    for p in packets:
        total_cmf_energy += p.energy_cmf
        total_rf_energy += p.energy_rf

    print("Total CMF energy")
    print(total_cmf_energy)

    # Below is the Artis compensation for their method of packet rejection
    """
    energy_ratio = total_energy.sum().sum() / total_cmf_energy

    print("Energy ratio")
    print(energy_ratio)
    
    for p in packets:
        p.energy_cmf *= energy_ratio
        p.energy_rf *= energy_ratio

    for e in energy_df_rows:
        e *= energy_ratio
    
    for row in energy_plot_df_rows:
        row[1] *= energy_ratio
    """
    print("Total RF energy")
    print(total_rf_energy)

    energy_bins = np.logspace(2, 3.8, spectrum_bins)
    energy_out = np.zeros((len(energy_bins - 1), time_steps))

    # Process packets
    (
        energy_df_rows,
        energy_plot_df_rows,
        energy_out,
        deposition_estimator,
    ) = gamma_packet_loop(
        packets,
        grey_opacity,
        photoabsorption_opacity,
        pair_creation_opacity,
        electron_number_density_time,
        mass_density_time,
        inv_volume_time,
        iron_group_fraction_per_shell.to_numpy(),
        inner_velocities,
        outer_velocities,
        times,
        dt_array,
        effective_time_array,
        energy_bins,
        energy_df_rows,
        energy_plot_df_rows,
        energy_out,
    )

    # DataFrame of energy information
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

    # DataFrame of positron energies
    energy_plot_positrons = pd.DataFrame(
        data=energy_plot_positron_rows,
        columns=[
            "packet_index",
            "energy_input",
            "energy_input_r",
            "energy_input_time",
        ],
    )

    # DataFrame of estimated deposition
    # Multiply dataframes by inv_volume_time array
    # if per unit volume is needed
    energy_estimated_deposition = (
        pd.DataFrame(data=deposition_estimator, columns=times[:-1])
    ) / dt_array

    # Energy is eV/s
    energy_df = pd.DataFrame(data=energy_df_rows, columns=times[:-1]) / dt_array

    final_energy = 0
    for p in packets:
        final_energy += p.energy_rf

    print("Final energy to test for conservation")
    print(final_energy)

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
