import numpy as np
from tqdm.auto import tqdm
import pandas as pd
import astropy.units as u
from numba import njit
from numba.typed import List

from tardis.montecarlo.montecarlo_numba import njit_dict_no_parallel
from tardis.energy_input.GXPacket import GXPacket, GXPacketStatus
from tardis.transport.geometry.calculate_distances_3d import (
    distance_trace_gamma
)
from tardis.transport.r_packet_transport_3d import (
    move_packet_3d
)
from tardis.energy_input.energy_source import (
    ni56_chain_energy,
    read_artis_lines,
)
from tardis.energy_input.calculate_opacity import (
    compton_opacity_calculation,
    photoabsorption_opacity_calculation,
    pair_creation_opacity_calculation,
    photoabsorption_opacity_calculation_kasen,
    pair_creation_opacity_artis,
    SIGMA_T,
)
from tardis.energy_input.gamma_ray_interactions import (
    get_compton_fraction_artis,
    scatter_type,
    compton_scatter,
    pair_creation_packet,
)
from tardis.energy_input.util import (
    doppler_gamma,
    C_CGS,
    H_CGS_KEV,
    kappa_calculation,
    get_index,
    get_random_unit_vector,
)
from tardis import constants as const
from tardis.energy_input.gamma_ray_estimators import deposition_estimator_kasen

# Energy: keV, exported as eV for SF solver
# distance: cm
# mass: g
# time: s


@njit(**njit_dict_no_parallel)
def sample_energy(energy, intensity):
    """Samples energy from energy and intensity

    Parameters
    ----------
    energy :  One-dimensional Numpy Array, dtype float
        Array of energies
    intensity :  One-dimensional Numpy Array, dtype float
        Array of intensities

    Returns
    -------
    float
        Energy
    """
    z = np.random.random()

    average = (energy * intensity).sum()
    total = 0
    for (e, i) in zip(energy, intensity):
        total += e * i / average
        if z <= total:
            return e

    return False


@njit(**njit_dict_no_parallel)
def sample_decay_time(isotope, ni56_tau, co56_tau):
    """Samples the decay time from the mean half-life
    of the isotopes (needs restructuring for more isotopes)

    Parameters
    ----------
    isotope : string
        Isotope as e.g. Ni56
    taus : dict
        Mean half-life for each isotope

    Returns
    -------
    float
        Decay time in seconds
    """
    if isotope == "Ni56":
        decay_time = -ni56_tau * np.log(np.random.random())
    else:
        decay_time = -ni56_tau * np.log(np.random.random()) - co56_tau * np.log(
            np.random.random()
        )
    return decay_time


@njit(**njit_dict_no_parallel)
def initialize_packets(
    decays_per_shell,
    input_energy,
    ni56_lines,
    co56_lines,
    inner_velocities,
    outer_velocities,
    inv_volume_time,
    times,
    energy_df_rows,
    effective_times,
    ni56_tau,
    co56_tau,
):
    """Initialize a list of GXPacket objects for the simulation
    to operate on.

    Parameters
    ----------
    decays_per_shell : array int64
        Number of decays per simulation shell
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
    ni56_tau, co56_tau : float64
        Mean half-life for each isotope

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

    number_of_packets = np.sum(decays_per_shell)

    print("Total packets:", number_of_packets)

    packet_energy = input_energy / number_of_packets

    print("Energy per packet", packet_energy)

    ni56_energy = (ni56_lines[:, 0] * ni56_lines[:, 1]).sum()
    co56_energy = (co56_lines[:, 0] * co56_lines[:, 1]).sum()

    ni56_fraction = ni56_energy / (ni56_energy + co56_energy)

    energy_plot_df_rows = np.zeros((number_of_packets, 8))
    energy_plot_positron_rows = np.zeros((number_of_packets, 4))

    j = 0
    for k, shell in enumerate(decays_per_shell):
        z = np.random.random(shell)

        initial_radii = (
            z * inner_velocities[k] ** 3.0
            + (1.0 - z) * outer_velocities[k] ** 3.0
        ) ** (1.0 / 3.0)

        for i in range(shell):
            decay_time = np.inf
            while decay_time > times[-1]:
                ni_or_co_random = np.random.random()

                if ni_or_co_random > ni56_fraction:
                    decay_type = "Co56"
                    energy = co56_lines[:, 0] * 1000
                    intensity = co56_lines[:, 1]
                    # positron energy scaled by intensity
                    positron_energy = 0.63 * 1000 * 0.19
                    positron_fraction = positron_energy / np.sum(
                        energy * intensity
                    )
                else:
                    decay_type = "Ni56"
                    energy = ni56_lines[:, 0] * 1000
                    intensity = ni56_lines[:, 1]
                    positron_fraction = 0

                decay_time = sample_decay_time(decay_type, ni56_tau, co56_tau)

            cmf_energy = sample_energy(energy, intensity)

            energy_factor = 1.0
            if decay_time < times[0]:
                energy_factor = decay_time / times[0]
                decay_time = times[0]

            decay_time_index = get_index(decay_time, times)

            energy_df_rows[k, decay_time_index] += (
                positron_fraction * packet_energy * 1000
            )

            scaled_r = initial_radii[i] * effective_times[decay_time_index]

            initial_location = scaled_r * get_random_unit_vector()
            initial_direction = get_random_unit_vector()
            initial_energy = packet_energy * energy_factor

            # draw a random gamma-ray in shell
            packet = GXPacket(
                initial_location,
                initial_direction,
                1.0,
                initial_energy,
                0.0,
                0.0,
                GXPacketStatus.IN_PROCESS,
                k,
                decay_time,
            )

            packet.energy_rf = packet.energy_cmf / doppler_gamma(
                packet.direction,
                packet.location,
                packet.time_current,
            )

            energy_plot_df_rows[j] = np.array(
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

            energy_plot_positron_rows[j] = [
                j,
                positron_fraction * packet_energy * 1000,
                # * inv_volume_time[packet.shell, decay_time_index],
                packet.get_location_r(),
                packet.time_current,
            ]

            packet.nu_cmf = cmf_energy / H_CGS_KEV

            packet.nu_rf = packet.nu_cmf / doppler_gamma(
                packet.direction,
                packet.location,
                packet.time_current,
            )

            packets.append(packet)

            j += 1

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
    path_to_artis_lines="",
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
    raw_isotope_abundance = model.raw_isotope_abundance

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
        plasma.number_density.loc[
            atom_number
        ] = plasma.isotope_number_density.loc[atom_number].values

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

    # Calculate number of packets per shell based on the mass of Ni56
    mass_ni56 = raw_isotope_abundance.loc[(28, 56)] * shell_masses
    number_ni56 = mass_ni56 / 56 * const.N_A
    total_number_ni56 = number_ni56.sum()

    decayed_packet_count = np.zeros(len(number_ni56), dtype=np.int64)

    for i, shell in enumerate(number_ni56):
        decayed_packet_count[i] = round(num_decays * shell / total_number_ni56)

    # hardcoded Ni and Co 56 mean half-lives
    taus = {
        "Ni56": 6.075 * u.d.to("s") / np.log(2),
        "Co56": 77.233 * u.d.to("s") / np.log(2),
    }

    # This will use the decay radiation database and be a more complex network eventually
    ni56_lines = read_artis_lines("ni56", path_to_artis_lines)
    co56_lines = read_artis_lines("co56", path_to_artis_lines)

    # urilight chooses to have 0 as the baseline for this calculation
    # but time_start should also be valid
    total_energy = ni56_chain_energy(
        taus, 0, time_end, number_ni56, ni56_lines, co56_lines
    )

    print("Total gamma-ray energy")
    print(total_energy.sum().sum() * u.keV.to("erg"))

    print("Total positron energy")
    print(total_energy["Co56"].sum(axis=0) * 0.0337 * u.keV.to("erg"))

    # Taking iron group to be elements 21-30
    # Used as part of the approximations for photoabsorption and pair creation
    # Dependent on atomic data
    iron_group_fraction_per_shell = model.abundance.loc[(21):(30)].sum(axis=0)

    # Need to update volume for positron deposition to be time-dependent
    print("Initializing packets")
    (
        packets,
        energy_df_rows,
        energy_plot_df_rows,
        energy_plot_positron_rows,
    ) = initialize_packets(
        decayed_packet_count,
        total_energy.sum().sum(),
        ni56_lines.to_numpy(),
        co56_lines.to_numpy(),
        inner_velocities,
        outer_velocities,
        inv_volume_time,
        times,
        energy_df_rows,
        effective_time_array,
        taus["Ni56"],
        taus["Co56"],
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

    energy_bins = np.logspace(2, 4, spectrum_bins)
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


@njit(**njit_dict_no_parallel)
def process_packet_path(packet):
    """Move the packet through interactions

    Parameters
    ----------
    packet : GXPacket
        Packet for processing

    Returns
    -------
    GXPacket
        Packet after processing
    float
        Energy injected into the ejecta
    """
    if packet.status == GXPacketStatus.COMPTON_SCATTER:
        comoving_freq_energy = packet.nu_cmf * H_CGS_KEV

        compton_angle, compton_fraction = get_compton_fraction_artis(
            comoving_freq_energy
        )

        # Packet is no longer a gamma-ray, destroy it
        if np.random.random() < 1 / compton_fraction:
            packet.nu_cmf = packet.nu_cmf / compton_fraction

            packet.direction = compton_scatter(packet, compton_angle)

            # Calculate rest frame frequency after scaling by the fraction that remains
            doppler_factor = doppler_gamma(
                packet.direction,
                packet.location,
                packet.time_current,
            )

            packet.nu_rf = packet.nu_cmf / doppler_factor
            packet.energy_rf = packet.energy_cmf / doppler_factor

            ejecta_energy_gained = 0.0
        else:
            packet.status = GXPacketStatus.PHOTOABSORPTION

    if packet.status == GXPacketStatus.PAIR_CREATION:
        packet = pair_creation_packet(packet)
        ejecta_energy_gained = 0.0

    if packet.status == GXPacketStatus.PHOTOABSORPTION:
        # Ejecta gains comoving energy
        ejecta_energy_gained = packet.energy_cmf

    return packet, ejecta_energy_gained


@njit(**njit_dict_no_parallel)
def gamma_packet_loop(
    packets,
    grey_opacity,
    photoabsorption_opacity_type,
    pair_creation_opacity_type,
    electron_number_density_time,
    mass_density_time,
    inv_volume_time,
    iron_group_fraction_per_shell,
    inner_velocities,
    outer_velocities,
    times,
    dt_array,
    effective_time_array,
    energy_bins,
    energy_df_rows,
    energy_plot_df_rows,
    energy_out,
):
    """Propagates packets through the simulation

    Parameters
    ----------
    packets : list
        List of GXPacket objects
    grey_opacity : float
        Grey opacity value in cm^2/g
    electron_number_density_time : array float64
        Electron number densities with time
    mass_density_time : array float64
        Mass densities with time
    inv_volume_time : array float64
        Inverse volumes with time
    iron_group_fraction_per_shell : array float64
        Iron group fraction per shell
    inner_velocities : array float64
        Inner velocities of the shells
    outer_velocities : array float64
        Inner velocities of the shells
    times : array float64
        Simulation time steps
    dt_array : array float64
        Simulation delta-time steps
    effective_time_array : array float64
        Simulation middle time steps
    energy_bins : array float64
        Bins for escaping gamma-rays
    energy_df_rows : array float64
        Energy output
    energy_plot_df_rows : array float64
        Energy output for plotting
    energy_out : array float64
        Escaped energy array

    Returns
    -------
    array float64
        Energy output
    array float64
        Energy output for plotting
    array float64
        Escaped energy array

    Raises
    ------
    ValueError
        Packet time index less than zero
    """
    escaped_packets = 0
    scattered_packets = 0
    packet_count = len(packets)
    print("Entering gamma ray loop for " + str(packet_count) + " packets")

    deposition_estimator = np.zeros_like(energy_df_rows)

    for i in range(packet_count):
        packet = packets[i]
        time_index = get_index(packet.time_current, times)

        if time_index < 0:
            print(packet.time_current, time_index)
            raise ValueError("Packet time index less than 0!")

        scattered = False

        initial_energy = packet.energy_cmf

        while packet.status == GXPacketStatus.IN_PROCESS:
            # Get delta-time value for this step
            dt = dt_array[time_index]
            # Calculate packet comoving energy for opacities
            comoving_energy = H_CGS_KEV * packet.nu_cmf

            if grey_opacity < 0:
                doppler_factor = doppler_gamma(
                    packet.direction,
                    packet.location,
                    effective_time_array[time_index],
                )

                kappa = kappa_calculation(comoving_energy)

                # artis threshold for Thomson scattering
                if kappa < 1e-2:
                    compton_opacity = (
                        SIGMA_T
                        * electron_number_density_time[packet.shell, time_index]
                    )
                else:
                    compton_opacity = compton_opacity_calculation(
                        comoving_energy,
                        electron_number_density_time[packet.shell, time_index],
                    )

                if photoabsorption_opacity_type == "kasen":
                    # currently not functional, requires proton count and
                    # electron count per isotope
                    photoabsorption_opacity = 0
                    # photoabsorption_opacity_calculation_kasen()
                else:
                    photoabsorption_opacity = (
                        photoabsorption_opacity_calculation(
                            comoving_energy,
                            mass_density_time[packet.shell, time_index],
                            iron_group_fraction_per_shell[packet.shell],
                        )
                    )

                if pair_creation_opacity_type == "artis":
                    pair_creation_opacity = pair_creation_opacity_artis(
                        comoving_energy,
                        mass_density_time[packet.shell, time_index],
                        iron_group_fraction_per_shell[packet.shell],
                    )
                else:
                    pair_creation_opacity = pair_creation_opacity_calculation(
                        comoving_energy,
                        mass_density_time[packet.shell, time_index],
                        iron_group_fraction_per_shell[packet.shell],
                    )

            else:
                compton_opacity = 0.0
                pair_creation_opacity = 0.0
                photoabsorption_opacity = (
                    grey_opacity * mass_density_time[packet.shell, time_index]
                )

            # convert opacities to rest frame
            total_opacity = (
                compton_opacity
                + photoabsorption_opacity
                + pair_creation_opacity
            ) * doppler_factor

            packet.tau = -np.log(np.random.random())

            (
                distance_interaction,
                distance_boundary,
                distance_time,
                shell_change,
            ) = distance_trace_gamma(
                packet,
                inner_velocities,
                outer_velocities,
                total_opacity,
                effective_time_array[time_index],
                times[time_index + 1],
            )

            distance = min(
                distance_interaction, distance_boundary, distance_time
            )

            packet.time_current += distance / C_CGS

            packet = move_packet_3d(packet, distance)

            deposition_estimator[packet.shell, time_index] += (
                (initial_energy * 1000)
                * distance
                * (packet.energy_cmf / initial_energy)
                * deposition_estimator_kasen(
                    comoving_energy,
                    mass_density_time[packet.shell, time_index],
                    iron_group_fraction_per_shell[packet.shell],
                )
            )

            if distance == distance_time:
                time_index += 1

                if time_index > len(effective_time_array) - 1:
                    # Packet ran out of time
                    packet.status = GXPacketStatus.END
                else:
                    packet.shell = get_index(
                        packet.get_location_r(),
                        inner_velocities * effective_time_array[time_index],
                    )

            elif distance == distance_interaction:

                packet.status = scatter_type(
                    compton_opacity,
                    photoabsorption_opacity,
                    total_opacity,
                )

                packet, ejecta_energy_gained = process_packet_path(packet)

                # Save packets to dataframe rows
                # convert KeV to eV / s / cm^3
                energy_df_rows[packet.shell, time_index] += (
                    ejecta_energy_gained * 1000
                )

                energy_plot_df_rows[i] = np.array(
                    [
                        i,
                        ejecta_energy_gained * 1000
                        # * inv_volume_time[packet.shell, time_index]
                        / dt,
                        packet.get_location_r(),
                        packet.time_current,
                        packet.shell,
                        compton_opacity,
                        photoabsorption_opacity,
                        pair_creation_opacity,
                    ]
                )

                if packet.status == GXPacketStatus.PHOTOABSORPTION:
                    # Packet destroyed, go to the next packet
                    break
                else:
                    packet.status = GXPacketStatus.IN_PROCESS
                    scattered = True

            else:
                packet.shell += shell_change

                if packet.shell > len(mass_density_time[:, 0]) - 1:
                    rest_energy = packet.nu_rf * H_CGS_KEV
                    bin_index = get_index(rest_energy, energy_bins)
                    bin_width = (
                        energy_bins[bin_index + 1] - energy_bins[bin_index]
                    )
                    energy_out[bin_index, time_index] += rest_energy / (
                        bin_width * dt
                    )
                    packet.status = GXPacketStatus.END
                    escaped_packets += 1
                    if scattered:
                        scattered_packets += 1
                elif packet.shell < 0:
                    packet.energy_rf = 0.0
                    packet.energy_cmf = 0.0
                    packet.status = GXPacketStatus.END

    print("Escaped packets:", escaped_packets)
    print("Scattered packets:", scattered_packets)

    return energy_df_rows, energy_plot_df_rows, energy_out, deposition_estimator
