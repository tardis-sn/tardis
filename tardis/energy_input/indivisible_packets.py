from bdb import effective
from random import sample
from time import time
import numpy as np
from tqdm.auto import tqdm
import pandas as pd
import astropy.units as u
from numba import njit, prange
from numba.typed import List
import radioactivedecay as rd

from tardis.energy_input.GXPacket import GXPacket, GXPacketStatus
from tardis.energy_input.gamma_ray_grid import (
    distance_trace,
    move_packet,
    read_artis_lines,
    mass_fraction_packets_per_shell,
    activity_per_shell,
    get_decay_database,
    get_tau,
    get_isotope_string,
    get_isotope
)
from tardis.energy_input.energy_source import (
    decay_nuclides,
    sample_mass,
    ni56_chain_energy,
    ni56_chain_energy_choice
)
from tardis.energy_input.calculate_opacity import (
    compton_opacity_calculation,
    photoabsorption_opacity_calculation,
    photoabsorption_opacity_calculation_kasen,
    pair_creation_opacity_calculation,
    SIGMA_T,
)
from tardis.energy_input.gamma_ray_interactions import (
    scatter_type,
    compton_scatter,
    pair_creation_packet,
    get_compton_fraction,
)
from tardis.energy_input.util import (
    get_random_theta_photon_array,
    get_random_phi_photon_array,
    doppler_gamma,
    BOUNDARY_THRESHOLD,
    C_CGS,
    H_CGS_KEV,
    kappa_calculation,
    get_index
)
from tardis import constants as const

# Energy: keV, exported as eV for SF solver
# distance: cm
# mass: g
# time: s


@njit
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
    dtype float
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

def sample_decay_time(isotope, taus):
    if isotope == "Ni56":
        decay_time = -taus["Ni56"] * np.log(np.random.random())
    else:
        decay_time = -taus["Ni56"] * np.log(np.random.random()) - taus["Co56"]  * np.log(np.random.random())
    return decay_time

def initialize_packets_equal_time(
    packet_count, 
    total_energy,
    power_df,
    times,
    shell_masses, 
    inner_velocities, 
    outer_velocities,
    ni56_lines,
    co56_lines,
    energy_df_rows,
    effective_time,
    taus,
    ):
    base_packet_energy = total_energy.sum().sum() / packet_count
    packets = List()

    energy_plot_df_rows = np.zeros((packet_count, 9))
    energy_plot_positron_rows = np.zeros((packet_count, 4))

    theta_locations = get_random_theta_photon_array(n=packet_count)
    phi_locations = get_random_phi_photon_array(n=packet_count)

    theta_directions = get_random_theta_photon_array(n=packet_count)
    phi_directions = get_random_phi_photon_array(n=packet_count)

    for i in tqdm(range(packet_count)):
        z = np.random.random()
        decay_time = z * times[0] + (1-z) * times[-1]
        decay_time_index = get_index(decay_time, times)
        
        time_df = power_df.loc[effective_time[decay_time_index]]

        # Need to limit sample space to isotope region
        # Slice with where?
        packet_radius, packet_shell = sample_mass(
            shell_masses[time_df.iloc[:, 0] > 0], 
            inner_velocities[time_df.iloc[:, 0] > 0], 
            outer_velocities[time_df.iloc[:, 0] > 0]
            )
        
        # randomly determine if Ni56 or Co56 for this shell at this time
        packet_isotope = get_isotope(time_df, packet_shell)

        if packet_isotope == "Ni56":
            # factor of 1000 for MeV -> keV
            energy = ni56_lines.energy.to_numpy() * 1000
            intensity = ni56_lines.intensity.to_numpy()
            positron_fraction = 0

        if packet_isotope == "Co56":
            energy = co56_lines.energy.to_numpy() * 1000
            intensity = co56_lines.intensity.to_numpy()
            # positron energy scaled by intensity
            positron_energy = 0.63 * 1000 * 0.19
            positron_fraction = positron_energy / (energy * intensity).sum()

        # average energy over time
        average_energy = ni56_chain_energy_choice(
            taus, 
            times[0], 
            times[-1], 
            time_df[packet_isotope][packet_shell] / taus["Ni56"], 
            ni56_lines, 
            co56_lines,
            packet_isotope
            ) / (times[-1] - times[0])

        # energy in timestep
        current_energy = ni56_chain_energy_choice(
            taus, 
            times[0], 
            decay_time, 
            time_df[packet_isotope][packet_shell] / taus["Ni56"], 
            ni56_lines, 
            co56_lines,
            packet_isotope
            ) / (decay_time - times[0])

        # scale energy by ratio of current time step energy to average energy over time
        packet_energy = base_packet_energy * current_energy / average_energy

        cmf_energy = sample_energy(energy, intensity)

        packet = GXPacket(
            location_r=packet_radius * effective_time[decay_time_index],
            location_theta=theta_locations[i],
            location_phi=phi_locations[i],
            direction_theta=theta_directions[i],
            direction_phi=phi_directions[i],
            energy_rf=1,
            energy_cmf=packet_energy,
            status=GXPacketStatus.IN_PROCESS,
            shell=packet_shell,
            nu_rf=0,
            nu_cmf=0,
            time_current=decay_time
        )

        packet.energy_rf = packet.energy_cmf / doppler_gamma(
                    packet.get_direction_vector(),
                    packet.get_position_vector(),
                    packet.time_current
                )

        packet.nu_cmf = cmf_energy / H_CGS_KEV

        packet.nu_rf = packet.nu_cmf / doppler_gamma(
            packet.get_direction_vector(),
            packet.get_position_vector(),
            packet.time_current
        )

        packets.append(packet)

        energy_df_rows[packet_shell, decay_time_index] += (
            positron_fraction
            * packet_energy
            * 1000
            #/ shell_volume
        )

        energy_plot_df_rows[i] = np.array(
            [
                i,
                packet.energy_cmf,
                packet.location_r,
                packet.location_theta,
                packet.time_current,
                int(packet.status),
                0,
                0,
                0,
            ]
        )

        energy_plot_positron_rows[i] = [
            i,
            positron_fraction
            * packet_energy
            * 1000,
            #/ shell_volume,
            packet.location_r,
            packet.time_current,
        ]


    return packets, energy_df_rows, energy_plot_df_rows, energy_plot_positron_rows




def initialize_packets(
    number_of_shells,
    decays_per_shell,
    model,
    inner_velocities,
    outer_velocities,
    ni56_lines,
    co56_lines,
    taus,
    input_energy_df,
    energy_df_rows,
    times,
    effective_times
):
    """Initializes packet properties
    and appends beta decay energy to output tables

    Parameters
    ----------
    number_of_shells : int
        Number of shells in model
    decays_per_shell : pandas.DataFrame
        Number of decays in a shell
    ejecta_volume : numpy.array
        Volume per shell
    inner_velocities : numpy.array
        Shell inner velocities
    outer_velocities : numpy.array
        Shell outer velocities
    decay_rad_db : pandas.DataFrame
        Decay radiation database
    scaled_activity_df : pandas.DataFrame
        Activity scaled per shell per isotope

    Returns
    -------
    list
        GXpacket objects
    numpy array
        energy binned per shell
    list
        packet info
    """
    packets = List()

    packet_energy = input_energy_df.sum().sum() / decays_per_shell.sum().sum()

    energy_plot_df_rows = np.zeros((decays_per_shell.sum(axis=0).sum(), 9))
    energy_plot_positron_rows = np.zeros((decays_per_shell.sum(axis=0).sum(), 4))

    j = 0
    for column in decays_per_shell:

        energy_per_shell = input_energy_df[column].to_numpy()
        packets_per_shell = decays_per_shell[column].to_numpy()

        if column == "Ni56":
            # factor of 1000 for MeV -> keV
            energy = ni56_lines.energy.to_numpy() * 1000
            intensity = ni56_lines.intensity.to_numpy()
            positron_fraction = 0

        if column == "Co56":
            energy = co56_lines.energy.to_numpy() * 1000
            intensity = co56_lines.intensity.to_numpy()
            # positron energy scaled by intensity
            positron_energy = 0.63 * 1000 * 0.19
            positron_fraction = positron_energy / (energy * intensity).sum()

        total_packets_for_isotope = packets_per_shell.sum()
        cmf_energy_array = np.zeros(total_packets_for_isotope)

        for i, _ in enumerate(cmf_energy_array):
            cmf_energy_array[i] = sample_energy(energy, intensity)

            if not cmf_energy_array[i]:
                print("No energy selected for this gamma ray!")
                continue

        for shell in range(number_of_shells):
            shell_energy = energy_per_shell[shell]
            shell_packets = packets_per_shell[shell]

            z = np.random.random(shell_packets)

            initial_radii = (
                z * inner_velocities[shell] ** 3.0
                + (1.0 - z) * outer_velocities[shell] ** 3.0
            ) ** (1.0 / 3.0)

            theta_locations = get_random_theta_photon_array(n=shell_packets)
            phi_locations = get_random_phi_photon_array(n=shell_packets)

            theta_directions = get_random_theta_photon_array(n=shell_packets)
            phi_directions = get_random_phi_photon_array(n=shell_packets)

            for i in range(shell_packets):
                while True:
                    decay_time = sample_decay_time(column, taus)
                    if decay_time < times[-1]:
                        break

                energy_factor = 1.0
                if decay_time < times[0]:
                    energy_factor = decay_time/times[0]
                    decay_time = times[0]

                decay_time_index = np.searchsorted(times, decay_time, side="right") - 1
                #model.time_explosion = decay_time * u.s
                #shell_volume = model.volume[shell].value

                # Add positron energy to the medium
                # convert KeV to eV / cm^3
                energy_df_rows[shell, decay_time_index] += (
                    positron_fraction
                    * packet_energy
                    * 1000
                    #/ shell_volume
                )

                # draw a random gamma-ray in shell
                packet = GXPacket(
                    location_r=initial_radii[i] * effective_times[decay_time_index],
                    location_theta=theta_locations[i],
                    location_phi=phi_locations[i],
                    direction_theta=theta_directions[i],
                    direction_phi=phi_directions[i],
                    energy_rf=1,
                    energy_cmf=packet_energy * energy_factor,
                    status=GXPacketStatus.IN_PROCESS,
                    shell=shell,
                    nu_rf=0,
                    nu_cmf=0,
                    time_current=decay_time
                )

                energy_plot_df_rows[j] = np.array(
                    [
                        i,
                        packet.energy_cmf,
                        packet.location_r,
                        packet.location_theta,
                        packet.time_current,
                        int(packet.status),
                        0,
                        0,
                        0,
                    ]
                )

                packet.energy_rf = packet.energy_cmf / doppler_gamma(
                    packet.get_direction_vector(),
                    packet.get_position_vector(),
                    packet.time_current
                )

                energy_plot_positron_rows[j] = [
                    j,
                    positron_fraction
                    * packet_energy
                    * 1000,
                    #/ shell_volume,
                    packet.location_r,
                    packet.time_current,
                ]

                packet.nu_cmf = cmf_energy_array[i] / H_CGS_KEV

                packet.nu_rf = packet.nu_cmf / doppler_gamma(
                    packet.get_direction_vector(),
                    packet.get_position_vector(),
                    packet.time_current
                )

                packets.append(packet)

                j += 1

    return packets, energy_df_rows, energy_plot_df_rows, energy_plot_positron_rows


def urilight_initialize_packets(
    decays_per_shell,
    input_energy_df,
    ni56_lines,
    co56_lines,
    inner_velocities,
    outer_velocities,
    times,
    energy_df_rows,
    effective_times,
    taus
):

    packets = List()

    number_of_packets = np.sum(decays_per_shell)

    print("Total packets:", number_of_packets)

    packet_energy = input_energy_df.sum().sum() / number_of_packets

    print("Energy per packet", packet_energy)

    ni56_energy = (ni56_lines.energy * ni56_lines.intensity).sum() 
    co56_energy = (co56_lines.energy * co56_lines.intensity).sum()

    ni56_fraction = ni56_energy / (ni56_energy + co56_energy)

    energy_plot_df_rows = np.zeros((number_of_packets, 9))
    energy_plot_positron_rows = np.zeros((number_of_packets, 4))

    theta_locations = get_random_theta_photon_array(n=number_of_packets)
    phi_locations = get_random_phi_photon_array(n=number_of_packets)

    theta_directions = get_random_theta_photon_array(n=number_of_packets)
    phi_directions = get_random_phi_photon_array(n=number_of_packets)

    j = 0
    for k, shell in tqdm(enumerate(decays_per_shell)):
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
                    energy = co56_lines.energy.to_numpy() * 1000
                    intensity = co56_lines.intensity.to_numpy()
                    # positron energy scaled by intensity
                    positron_energy = 0.63 * 1000 * 0.19
                    positron_fraction = positron_energy / (energy * intensity).sum()
                else:
                    decay_type = "Ni56"
                    energy = ni56_lines.energy.to_numpy() * 1000
                    intensity = ni56_lines.intensity.to_numpy()
                    positron_fraction = 0
                
                decay_time = sample_decay_time(decay_type, taus)

            cmf_energy = sample_energy(energy, intensity)

            energy_factor = 1.0
            if decay_time < times[0]:
                energy_factor = decay_time/times[0]
                decay_time = times[0]

            decay_time_index = get_index(decay_time, times)

            energy_df_rows[k, decay_time_index] += (
                positron_fraction
                * packet_energy
                * 1000
                #/ shell_volume
            )

            # draw a random gamma-ray in shell
            packet = GXPacket(
                location_r=initial_radii[i] * effective_times[decay_time_index],
                location_theta=theta_locations[j],
                location_phi=phi_locations[j],
                direction_theta=theta_directions[j],
                direction_phi=phi_directions[j],
                energy_rf=1,
                energy_cmf=packet_energy * energy_factor,
                status=GXPacketStatus.IN_PROCESS,
                shell=k,
                nu_rf=0,
                nu_cmf=0,
                time_current=decay_time
            )

            packet.energy_rf = packet.energy_cmf / doppler_gamma(
                packet.get_direction_vector(),
                packet.get_position_vector(),
                packet.time_current
            )

            energy_plot_df_rows[j] = np.array(
                [
                    i,
                    packet.energy_rf,
                    packet.location_r,
                    packet.location_theta,
                    packet.time_current,
                    int(packet.status),
                    0,
                    0,
                    0,
                ]
            )

            energy_plot_positron_rows[j] = [
                j,
                positron_fraction
                * packet_energy
                * 1000,
                #/ shell_volume,
                packet.location_r,
                packet.time_current,
            ]

            packet.nu_cmf = cmf_energy / H_CGS_KEV

            packet.nu_rf = packet.nu_cmf / doppler_gamma(
                packet.get_direction_vector(),
                packet.get_position_vector(),
                packet.time_current
            )

            packets.append(packet)

            j += 1

    return packets, energy_df_rows, energy_plot_df_rows, energy_plot_positron_rows



def main_gamma_ray_loop(num_decays, model, plasma, time_steps=10, time_end=80.0, grey_opacity=0.0):
    """Main loop that determines the gamma ray propagation

    Parameters
    ----------
    num_decays : int
        Number of decays requested
    model : tardis.Radial1DModel
        The tardis model to calculate gamma ray propagation through
    plasma : tardis.plasma.BasePlasma
        The tardis plasma with calculated atomic number density
    grey_opacity : float
        Grey photoabsorption opacity for gamma-rays in cm^2 g^-1

    Returns
    -------
    pandas.DataFrame
        Energy per mass per shell in units of eV/s/cm^-3
    pandas.DataFrame
        Columns:
        packet index,
        Energy input per packet,
        radius of deposition,
        theta angle of deposition,
        time of deposition,
        type of deposition where:
            -1 = beta decay,
            0 = Compton scatter,
            1 = photoabsorption,
            2 = pair creation
    list
        Energy of escaping packets
    numpy.ndarray
        Scaled activity per shell
    pandas.DataFrame
        Energy injected into the model per shell
    """
    escape_energy = []

    # Enforce cgs
    outer_velocities = model.v_outer.to("cm/s").value
    inner_velocities = model.v_inner.to("cm/s").value
    ejecta_density = model.density.to("g/cm^3").value
    ejecta_volume = model.volume.to("cm^3").value
    ejecta_velocity_volume = 4 * np.pi / 3 * (outer_velocities ** 3.0 - inner_velocities ** 3.0)
    time_explosion = model.time_explosion.to("s").value
    number_of_shells = model.no_of_shells
    raw_isotope_abundance = model.raw_isotope_abundance

    shell_masses = ejecta_volume * ejecta_density

    time_start = time_explosion
    time_end *= u.d.to(u.s)

    assert time_start < time_end, "Error, simulation start time greater than end time!"

    times = np.zeros(time_steps + 1)

    # log time steps
    for i in range(time_steps + 1):
        times[i]=np.log(time_start) + (np.log(time_end) - np.log(time_start)) / time_steps * i
        times[i]=np.exp(times[i])

    dt_array = np.diff(times)
    effective_time_array = np.array([np.sqrt(times[i] * times[i + 1]) for i in range(time_steps)])

    electron_number_density = (
        plasma.number_density.mul(plasma.number_density.index, axis=0)
    ).sum()

    electron_number_density_time = np.zeros((len(ejecta_velocity_volume), len(effective_time_array)))
    mass_density_time = np.zeros((len(ejecta_velocity_volume), len(effective_time_array)))
    electron_number = (electron_number_density * ejecta_volume).to_numpy()

    for i, t in enumerate(effective_time_array):
        mass_density_time[:, i] = shell_masses * (1.0 / ejecta_velocity_volume) / (t ** 3.0)
        electron_number_density_time[:, i] = electron_number * (1.0 / ejecta_velocity_volume) / (t ** 3.0)

    energy_df_rows = np.zeros((number_of_shells, time_steps))

    decay_count_rows = []

    """
    for time in effective_time_array:
        to_be_decayed = raw_isotope_abundance.copy()
        abundances_decayed = to_be_decayed.decay(time * u.s.to('d')).T
        dataframe_columns = [f"{rd.utils.Z_DICT[i[0]]}{i[1]}" for i in abundances_decayed.columns]
        abundances_decayed.columns = dataframe_columns

        decay_count_rows.append(mass_fraction_packets_per_shell(abundances_decayed.multiply(shell_masses, axis='index'), num_decays))
        #activity_df, decay, meta = activity_per_shell(abundances_decayed.multiply(shell_masses, axis='index'))
        #decay_count_rows.append(activity_df)
    """
    mass_ni56 = raw_isotope_abundance.loc[(28, 56)] * shell_masses
    number_ni56 = mass_ni56 / 56 * const.N_A
    total_number_ni56 = number_ni56.sum()

    decayed_packet_count = np.zeros(len(number_ni56), dtype=np.int64)

    for i, shell in enumerate(number_ni56):
        decayed_packet_count[i] = round(num_decays * shell / total_number_ni56)

    #decayed_packet_count = pd.concat(decay_count_rows, keys=effective_time_array, names=["time", "shell"])

    #decay_db, meta = get_decay_database(decayed_packet_count)

    taus = {"Ni56": 6.075 * u.d.to('s') / np.log(2), "Co56": 77.233 * u.d.to('s') / np.log(2)}

    """
    for column in decayed_packet_count.columns:
        if column == "Fe56":
            continue
        taus[column] = get_tau(meta, column).value
    """
    # This will use the decay radiation database and be a more complex network eventually

    ni56_lines = read_artis_lines("ni56")
    co56_lines = read_artis_lines("co56")

    # urilight chooses to have 0 as the baseline for this calculation
    # but time_start should also be valid
    total_energy = ni56_chain_energy(taus, 0, time_end, number_ni56, ni56_lines, co56_lines)

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
    """
    packets, energy_df_rows, energy_plot_df_rows, energy_plot_positron_rows = initialize_packets(
        number_of_shells,
        decayed_packet_count,
        model,
        inner_velocities,
        outer_velocities,
        ni56_lines,
        co56_lines,
        taus,
        total_energy,
        energy_df_rows,
        times,
        effective_time_array
    )
    packets, energy_df_rows, energy_plot_df_rows, energy_plot_positron_rows = initialize_packets_equal_time(
        num_decays, 
        total_energy, 
        decayed_packet_count, 
        times, 
        shell_masses, 
        inner_velocities, 
        outer_velocities, 
        ni56_lines, 
        co56_lines, 
        energy_df_rows,
        effective_time_array,
        taus
        )
    """
    
    (packets, 
    energy_df_rows, 
    energy_plot_df_rows, 
    energy_plot_positron_rows) = urilight_initialize_packets(
        decayed_packet_count, 
        total_energy, 
        ni56_lines, 
        co56_lines, 
        inner_velocities, 
        outer_velocities, 
        times, 
        energy_df_rows, 
        effective_time_array, 
        taus
        )

    print("Total positron energy from packets")
    print((energy_df_rows).sum().sum() * u.eV.to("erg"))

    #energy_plot_positron_rows[:, 1] /= dt_array[0]

    total_cmf_energy = 0
    total_rf_energy = 0

    for p in packets:
        total_cmf_energy += p.energy_cmf
        total_rf_energy += p.energy_rf

    energy_ratio = total_energy.sum().sum() / total_cmf_energy

    print("Total CMF energy")
    print(total_cmf_energy)

    print("Energy ratio")
    print(energy_ratio)
    """
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

    # Need to update volume with time-step for deposition to be time-dependent
    energy_df_rows, energy_plot_df_rows = gamma_packet_loop(
        packets,
        grey_opacity,
        electron_number_density_time,
        mass_density_time,
        iron_group_fraction_per_shell.to_numpy(),
        inner_velocities,
        outer_velocities,
        times,
        dt_array,
        effective_time_array,
        energy_df_rows,
        energy_plot_df_rows,
    )

    # DataFrame of energy information
    energy_plot_df = pd.DataFrame(
        data=energy_plot_df_rows,
        columns=[
            "packet_index",
            "energy_input",
            "energy_input_r",
            "energy_input_theta",
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
            "energy_input_time"
        ],
    )

    # Energy is eV/s/cm^-3
    energy_df = pd.DataFrame(data=energy_df_rows, columns=times[:-1]) / dt_array

    final_energy = 0
    for p in packets:
        final_energy += p.energy_rf

    print("Final energy to test for conservation")
    print(final_energy)

    return (
        energy_df,
        energy_plot_df,
        escape_energy,
        decayed_packet_count,
        energy_plot_positrons,
    )


@njit
def process_packet_path(packet):

    if packet.status == GXPacketStatus.COMPTON_SCATTER:
        comoving_freq_energy = packet.nu_cmf * H_CGS_KEV

        compton_angle, compton_fraction = get_compton_fraction(
            comoving_freq_energy
        )

        kappa = kappa_calculation(comoving_freq_energy)

        # Basic check to see if Thomson scattering needs to be considered later
        if kappa < 1e-2:
            print("Thomson should happen")

        # Packet is no longer a gamma-ray, destroy it
        if np.random.random() < compton_fraction:
            packet.status = GXPacketStatus.PHOTOABSORPTION
        else:
            ejecta_energy_gained = 0.0

        packet.nu_cmf = packet.nu_cmf * compton_fraction

        (
            packet.direction_theta,
            packet.direction_phi,
        ) = compton_scatter(packet, compton_angle)

        # Calculate rest frame frequency after scaling by the fraction that remains
        doppler_factor = doppler_gamma(
            packet.get_direction_vector(), packet.get_position_vector(), packet.time_current
        )

        packet.nu_rf = packet.nu_cmf / doppler_factor
        packet.energy_rf = packet.energy_cmf / doppler_factor

    if packet.status == GXPacketStatus.PAIR_CREATION:
        packet = pair_creation_packet(packet)
        ejecta_energy_gained = 0.0

    if packet.status == GXPacketStatus.PHOTOABSORPTION:
        # Ejecta gains comoving energy
        ejecta_energy_gained = packet.energy_cmf

    return packet, ejecta_energy_gained


@njit(parallel=True)
def gamma_packet_loop(
    packets,
    grey_opacity,
    electron_number_density_time,
    mass_density_time,
    iron_group_fraction_per_shell,
    inner_velocities,
    outer_velocities,
    times,
    dt_array,
    effective_time_array,
    energy_df_rows,
    energy_plot_df_rows,
):
    escaped_packets = 0
    scattered_packets = 0
    packet_count = len(packets)
    print("Entering gamma ray loop for " + str(packet_count) + " packets")

    for i in prange(packet_count):
        packet = packets[i]
        time_index = get_index(packet.time_current, times)
        
        if time_index < 0:
            print(packet.time_current, time_index)
            raise ValueError("Packet time index less than 0!")

        scattered = False

        while packet.status == GXPacketStatus.IN_PROCESS:
            # Get delta-time value for this step
            dt = dt_array[time_index]

            # Calculate packet comoving energy for opacities
            comoving_energy = H_CGS_KEV * packet.nu_cmf

            doppler_factor = doppler_gamma(
                packet.get_direction_vector(), 
                packet.get_position_vector(), 
                effective_time_array[time_index]
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
                        electron_number_density_time[packet.shell, time_index]
                    )
            """
            photoabsorption_opacity = photoabsorption_opacity_calculation(
                    comoving_energy,
                    mass_density_time[packet.shell, time_index],
                    iron_group_fraction_per_shell[packet.shell],
                )
            """
            photoabsorption_opacity = photoabsorption_opacity_calculation_kasen(
                    comoving_energy,
                    electron_number_density_time[packet.shell, time_index]
                )  

            pair_creation_opacity =  pair_creation_opacity_calculation(
                    comoving_energy,
                    electron_number_density_time[packet.shell, time_index],
                    iron_group_fraction_per_shell[packet.shell],
                )

            if grey_opacity > 0:
                compton_opacity = 0.0
                pair_creation_opacity = 0.0
                photoabsorption_opacity = grey_opacity * electron_number_density_time[packet.shell, time_index]
            
            # convert opacities to rest frame
            total_opacity = (
                compton_opacity
                + photoabsorption_opacity
                + pair_creation_opacity
            ) * doppler_factor 

            packet.tau = -np.log(np.random.random())
            
            (distance_interaction, distance_boundary, distance_time) = distance_trace(
                packet,
                inner_velocities,
                outer_velocities,
                total_opacity,
                effective_time_array[time_index],
                times[time_index + 1]
            )

            distance = min(distance_interaction, distance_boundary, distance_time)

            packet.time_current += (
                    distance / C_CGS
                )

            packet = move_packet(packet, distance)

            if distance == distance_time:
                time_index += 1
                
                if time_index > len(effective_time_array) - 1:
                    # Packet ran out of time
                    packet.status = GXPacketStatus.END
                else:
                    packet.shell = get_index(packet.location_r, inner_velocities * effective_time_array[time_index])

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
                    #/ ejecta_volume[packet.shell] 
                )
                
                energy_plot_df_rows[i] = np.array(
                    [
                        i,
                        ejecta_energy_gained
                        * 1000
                        #/ ejecta_volume[packet.shell] 
                        / dt,
                        packet.location_r,
                        packet.location_theta,
                        packet.time_current,
                        int(packet.status),
                        compton_opacity,
                        photoabsorption_opacity,
                        total_opacity,
                    ]
                )
                
                if packet.status == GXPacketStatus.PHOTOABSORPTION:
                    # Packet destroyed, go to the next packet
                    break
                else:
                    packet.status = GXPacketStatus.IN_PROCESS
                    scattered = True

            else:
                # packet.tau -= total_opacity * distance_boundary * time_explosion
                # overshoot so that the gamma-ray is comfortably in the next shell
                #packet = move_packet(
                #    packet, distance# * (1 + BOUNDARY_THRESHOLD)
                #)
            
                packet.shell += 1
                
                #= np.searchsorted(
                #        outer_velocities * effective_time_array[time_index], 
                #        packet.location_r,
                #        side="left"
                #    ) - 1

                if packet.shell > len(mass_density_time[:, 0]) - 1:
                    # escape_energy.append(packet.energy_rf)
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

    return energy_df_rows, energy_plot_df_rows
