from random import sample
import numpy as np
from tqdm.auto import tqdm
import pandas as pd
import astropy.units as u
from numba import njit, prange
from numba.typed import List

from tardis.energy_input.GXPacket import GXPacket, GXPacketStatus
from tardis.energy_input.gamma_ray_grid import (
    distance_trace,
    move_packet,
    read_artis_lines,
    mass_fraction_packets_per_shell,
    get_decay_database,
    get_tau,
    get_isotope_string
)
from tardis.energy_input.energy_source import (
    decay_nuclides,
)
from tardis.energy_input.calculate_opacity import (
    compton_opacity_calculation,
    photoabsorption_opacity_calculation,
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
    times
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

    energy_plot_df_rows = np.zeros((decays_per_shell.sum(axis=0).sum(), 9))
    energy_plot_positron_rows = np.zeros((decays_per_shell.sum(axis=0).sum(), 4))

    j = 0
    for column in decays_per_shell:

        tau = taus[column]

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
        cmf_energy = np.zeros(total_packets_for_isotope)

        for i, _ in enumerate(cmf_energy):
            cmf_energy[i] = sample_energy(energy, intensity)

            if not cmf_energy[i]:
                print("No energy selected for this gamma ray!")
                continue

        for shell in range(number_of_shells):
            shell_energy = energy_per_shell[shell]
            shell_packets = packets_per_shell[shell]
            energy_cmf = 0

            if shell_packets > 0:
                energy_cmf = shell_energy / shell_packets

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
                if column == "Ni56":
                    decay_time = -tau * np.log10(np.random.random())
                else:
                    decay_time = -taus["Ni56"] * np.log10(np.random.random()) - tau * np.log10(np.random.random())
                decay_time_index = np.argmin(np.abs(times-decay_time))

                energy_factor = 1.0
                if decay_time < times[0]:
                    energy_factor = decay_time/times[0]
                    decay_time = times[0]

                model.time_explosion = decay_time * u.s
                shell_volume = model.volume[shell].value

                # Add positron energy to the medium
                # convert KeV to eV / cm^3
                energy_df_rows[shell, decay_time_index] += (
                    positron_fraction
                    * shell_energy
                    * 1000
                    / shell_volume
                )
                # draw a random gamma-ray in shell
                packet = GXPacket(
                    location_r=initial_radii[i],
                    location_theta=theta_locations[i],
                    location_phi=phi_locations[i],
                    direction_theta=theta_directions[i],
                    direction_phi=phi_directions[i],
                    energy_rf=1,
                    energy_cmf=energy_cmf * energy_factor,
                    status=GXPacketStatus.IN_PROCESS,
                    shell=shell,
                    nu_rf=0,
                    nu_cmf=0,
                    time_current=decay_time
                )

                packet.energy_rf = packet.energy_cmf / doppler_gamma(
                    packet.get_direction_vector(),
                    packet.location_r,
                )

                energy_plot_positron_rows[j] = [
                    j,
                    positron_fraction
                    * shell_energy
                    * 1000
                    / shell_volume,
                    packet.location_r,
                    packet.location_theta,
                ]

                packet.nu_cmf = cmf_energy[i] / H_CGS_KEV

                packet.nu_rf = packet.nu_cmf / doppler_gamma(
                    packet.get_direction_vector(),
                    packet.location_r,
                )

                packets.append(packet)

                j += 1

    return packets, energy_df_rows, energy_plot_df_rows, energy_plot_positron_rows


def main_gamma_ray_loop(num_decays, model, plasma, time_steps=10, time_start=2.0, grey_opacity=0.0):
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

    # Note the use of velocity as the radial coordinate
    # Enforce cgs
    outer_velocities = model.v_outer.to("cm/s").value
    inner_velocities = model.v_inner.to("cm/s").value
    ejecta_density = model.density.to("g/cm^3").value
    ejecta_volume = model.volume.to("cm^3").value
    time_explosion = model.time_explosion.to("s").value
    number_of_shells = model.no_of_shells
    raw_isotope_abundance = model.raw_isotope_abundance

    shell_masses = ejecta_volume * ejecta_density

    time_end = time_explosion
    time_start *= u.day.to("s")

    times = np.zeros(time_steps + 1)

    # log time steps
    for i in range(time_steps + 1):
        times[i]=np.log(time_start)+(np.log(time_end)-np.log(time_start))/time_steps*i
        times[i]=np.exp(times[i])

    print(times)

    electron_number_density = (
        plasma.number_density.mul(plasma.number_density.index, axis=0)
    ).sum()

    energy_df_rows = np.zeros((number_of_shells, time_steps + 1))

    new_abundance_rows = []

    for index, mass in enumerate(shell_masses):
        isotope_abundance = raw_isotope_abundance[index]
        isotope_dict = {}
        for (
            atom_number,
            atom_mass,
        ), abundance in isotope_abundance.iteritems():
            isotope_string = get_isotope_string(atom_number, atom_mass)
            isotope_dict[isotope_string] = abundance
        new_abundance_rows.append(
            decay_nuclides(
                mass * u.g.to("M_sun"),
                isotope_dict,
                np.array([time_start * u.s.to("day")]),
            )
        )

    decayed_isotope_abundance = pd.concat(new_abundance_rows).T
    decayed_isotope_abundance.columns = raw_isotope_abundance.columns

    # scale abundances, needs to be done per-isotope in future
    decayed_isotope_abundance *= raw_isotope_abundance.values
    decayed_isotope_abundance = decayed_isotope_abundance.T

    decays_per_shell = mass_fraction_packets_per_shell(shell_masses, decayed_isotope_abundance, num_decays)

    decay_db, meta = get_decay_database(decayed_isotope_abundance)

    taus = {}

    for column in decayed_isotope_abundance.columns:
        if column == "Fe56":
            continue
        taus[column] = get_tau(meta, column).value

    # This will use the decay radiation database and be a more complex network eventually

    ni56_lines = read_artis_lines("ni56")
    co56_lines = read_artis_lines("co56")

    number_ni56 = raw_isotope_abundance.loc[(28, 56)] * shell_masses / 56 * const.N_A

    # urilight chooses to have 0 as the baseline for this calculation
    # but time_start should also be valid
    total_ni56 = -taus["Ni56"] * (np.exp(-time_end / taus["Ni56"]) - np.exp(-0 / taus["Ni56"]))
    total_co56 = -taus["Co56"] * (np.exp(-time_end / taus["Co56"]) - np.exp(-0 / taus["Co56"]))

    total_energy = pd.DataFrame()
    
    total_energy["Ni56"] = number_ni56 * (
            (ni56_lines.energy * 1000 * ni56_lines.intensity).sum()
            / taus["Ni56"]
            *
            total_ni56
    )

    total_energy["Co56"] = number_ni56 * (
            (co56_lines.energy * 1000 * co56_lines.intensity).sum()
            / (taus["Ni56"] - taus["Co56"])
            * (total_ni56 - total_co56)
    )

    print("Total expected energy")
    print(total_energy.sum().sum() * u.keV.to("erg"))

    # Taking iron group to be elements 21-30
    # Used as part of the approximations for photoabsorption and pair creation
    # Dependent on atomic data
    iron_group_fraction_per_shell = model.abundance.loc[(21):(30)].sum(axis=0)

    # Need to update volume for positron deposition to be time-dependent
    print("Initializing packets")
    packets, energy_df_rows, energy_plot_df_rows, energy_plot_positron_rows = initialize_packets(
        number_of_shells,
        decays_per_shell,
        model,
        inner_velocities,
        outer_velocities,
        ni56_lines,
        co56_lines,
        taus,
        total_energy,
        energy_df_rows,
        times
    )

    dt = times[1] - times[0]

    energy_df_rows /= dt
    energy_plot_positron_rows[:, 1] /= dt

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

    for p in packets:
        p.energy_cmf *= energy_ratio
        p.energy_rf *= energy_ratio

    for e in energy_df_rows:
        e *= energy_ratio

    # for row in energy_plot_df_rows:
    #    row[1] *= energy_ratio

    print("Total RF energy")
    print(total_rf_energy)

    model.time_explosion = time_end

    # Need to update volume with time-step for deposition to be time-dependent
    energy_df_rows, energy_plot_df_rows = gamma_packet_loop(
        packets,
        grey_opacity,
        electron_number_density.to_numpy(),
        ejecta_density,
        ejecta_volume,
        iron_group_fraction_per_shell.to_numpy(),
        inner_velocities,
        outer_velocities,
        times,
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
            "energy_input_theta"
        ],
    )

    # Energy is eV/s/cm^-3
    energy_df = pd.DataFrame(data=energy_df_rows, columns=times)

    final_energy = 0
    for p in packets:
        final_energy += p.energy_rf

    print("Final energy to test for conservation")
    print(final_energy)

    return (
        energy_df,
        energy_plot_df,
        escape_energy,
        decays_per_shell,
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
            packet.get_direction_vector(), packet.location_r
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
    electron_number_density,
    ejecta_density,
    ejecta_volume,
    iron_group_fraction_per_shell,
    inner_velocities,
    outer_velocities,
    times,
    energy_df_rows,
    energy_plot_df_rows,
):
    packet_count = len(packets)
    print("Entering gamma ray loop for " + str(packet_count) + " packets")

    for i in prange(packet_count):
        packet = packets[i]

        time_index = np.argmin(np.abs(times-packet.time_current))

        dt = 1

        while packet.status == GXPacketStatus.IN_PROCESS:

            time = times[time_index]
            
            # Get delta-time value for this step
            if time_index < len(times) - 1:
                dt = times[time_index + 1] - time

            # Calculate packet comoving energy for opacities
            comoving_energy = H_CGS_KEV * packet.nu_cmf

            doppler_factor = doppler_gamma(
                packet.get_direction_vector(), packet.location_r
            )

            kappa = kappa_calculation(comoving_energy)

            # artis threshold for Thomson scattering
            if kappa < 1e-2:
                compton_opacity = (
                    SIGMA_T
                    * electron_number_density[packet.shell]
                    * doppler_factor
                )
            else:
                compton_opacity = (
                    compton_opacity_calculation(
                        comoving_energy, electron_number_density[packet.shell]
                    )
                    * doppler_factor
                )

            photoabsorption_opacity = (
                photoabsorption_opacity_calculation(
                    comoving_energy,
                    ejecta_density[packet.shell],
                    iron_group_fraction_per_shell[packet.shell],
                )
                * doppler_factor
            )

            pair_creation_opacity = (
                pair_creation_opacity_calculation(
                    comoving_energy,
                    ejecta_density[packet.shell],
                    iron_group_fraction_per_shell[packet.shell],
                )
                * doppler_factor
            )

            if grey_opacity > 0:
                compton_opacity = 0.0
                pair_creation_opacity = 0.0
                photoabsorption_opacity = (
                    grey_opacity * ejecta_density[packet.shell]
                )

            # convert opacities to rest frame
            total_opacity = (
                compton_opacity
                + photoabsorption_opacity
                + pair_creation_opacity
            )

            packet.tau = -np.log(np.random.random())

            (distance_interaction, distance_boundary, distance_time) = distance_trace(
                packet,
                inner_velocities,
                outer_velocities,
                total_opacity,
                times[time_index + 1],
            )

            distance = min(distance_interaction, distance_boundary, distance_time)

            packet.time_current += (
                    distance / C_CGS
                )

            if distance == distance_time:
                packet = move_packet(packet, distance)
                time_index += 1
                
            elif distance == distance_interaction:

                packet.status = scatter_type(
                    compton_opacity,
                    photoabsorption_opacity,
                    total_opacity,
                )

                packet = move_packet(packet, distance)
                packet.shell = np.searchsorted(
                    outer_velocities, packet.location_r, side="left"
                )

                packet, ejecta_energy_gained = process_packet_path(packet)

                # Save packets to dataframe rows
                # convert KeV to eV / s / cm^3
                energy_df_rows[packet.shell, time_index] += (
                    ejecta_energy_gained * 1000 / ejecta_volume[packet.shell] / dt
                )

                energy_plot_df_rows[i] = np.array(
                    [
                        i,
                        ejecta_energy_gained
                        * 1000
                        / ejecta_volume[packet.shell] / dt,
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

            else:
                # packet.tau -= total_opacity * distance_boundary * time_explosion
                # overshoot so that the gamma-ray is comfortably in the next shell
                packet = move_packet(
                    packet, distance * (1 + BOUNDARY_THRESHOLD)
                )
                packet.shell = np.searchsorted(
                    outer_velocities, packet.location_r, side="left"
                )

            if packet.shell > len(ejecta_density) - 1:
                # escape_energy.append(packet.energy_rf)
                packet.status = GXPacketStatus.END
            elif packet.shell < 0:
                packet.energy_rf = 0.0
                packet.energy_cmf = 0.0
                packet.status = GXPacketStatus.END
            elif time_index > len(times) - 2:
                # Packet ran out of time
                packet.status = GXPacketStatus.END

    return energy_df_rows, energy_plot_df_rows
