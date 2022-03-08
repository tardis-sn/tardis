from random import sample
import numpy as np
import copy
from tqdm.auto import tqdm
import pandas as pd
import astropy.units as u
from numba import njit

from tardis.energy_input.GXPacket import GXPacket, GXPacketStatus
from tardis.energy_input.gamma_ray_grid import (
    distance_trace,
    move_packet,
    compute_required_photons_per_shell_artis,
    read_artis_lines,
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
    get_random_theta_photon,
    get_random_phi_photon,
    get_random_theta_photon_array,
    get_random_phi_photon_array,
    doppler_gamma,
    BOUNDARY_THRESHOLD,
    C_CGS,
    H_CGS_KEV,
    kappa_calculation,
)
from tardis import constants as const
from tardis.util.base import (
    atomic_number2element_symbol,
)

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
    ejecta_volume,
    inner_velocities,
    outer_velocities,
    ni56_lines,
    co56_lines,
    scaled_activity_df,
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
    packets = []
    energy_df_rows = np.zeros(number_of_shells)
    energy_plot_df_rows = []

    scaled_decays_per_shell = decays_per_shell.copy()

    for column in scaled_decays_per_shell:

        if column == "Ni56":
            # factor of 1000 for MeV -> keV
            energy = ni56_lines.energy.to_numpy() * 1000
            intensity = ni56_lines.intensity.to_numpy()
            positron_energy = 0
            mean_energy = (
                ni56_lines.energy.to_numpy()
                * 1000
                * ni56_lines.intensity.to_numpy()
            ).sum()

        if column == "Co56":
            energy = co56_lines.energy.to_numpy() * 1000
            intensity = co56_lines.intensity.to_numpy()
            # positron energy scaled by intensity
            positron_energy = 0.63 * 1000 * 0.19
            mean_energy = (
                co56_lines.energy.to_numpy()
                * 1000
                * co56_lines.intensity.to_numpy()
            ).sum()

        for shell in range(number_of_shells):
            activity = scaled_activity_df[column].iloc[shell]
            energy_cmf = mean_energy * activity
            packets_per_shell = scaled_decays_per_shell[column].iloc[shell]

            z = np.random.random(packets_per_shell)

            initial_radii = (
                z * inner_velocities[shell] ** 3.0
                + (1.0 - z) * outer_velocities[shell] ** 3.0
            ) ** (1.0 / 3.0)

            theta_locations = get_random_theta_photon_array(n=packets_per_shell)
            phi_locations = get_random_phi_photon_array(n=packets_per_shell)

            theta_directions = get_random_theta_photon_array(
                n=packets_per_shell
            )
            phi_directions = get_random_phi_photon_array(n=packets_per_shell)

            for i in range(packets_per_shell):
                # draw a random gamma-ray in shell
                packet = GXPacket(
                    location_r=initial_radii[i],
                    location_theta=theta_locations[i],
                    location_phi=phi_locations[i],
                    direction_theta=theta_directions[i],
                    direction_phi=phi_directions[i],
                    energy_rf=1,
                    energy_cmf=energy_cmf,
                    status=GXPacketStatus.IN_PROCESS,
                    shell=shell,
                    nu_rf=0,
                    nu_cmf=0,
                    activity=activity,
                )

                packet.energy_rf = packet.energy_cmf / doppler_gamma(
                    packet.get_direction_vector(),
                    packet.location_r,
                )

                # Add positron energy to the medium
                # convert KeV to eV / cm^3
                energy_df_rows[shell] += (
                    positron_energy * activity * 1000 / ejecta_volume[shell],
                )
                energy_plot_df_rows.append(
                    [
                        -1,
                        positron_energy
                        * activity
                        * 1000
                        / ejecta_volume[shell],
                        packet.location_r,
                        packet.location_theta,
                        0.0,
                        -1,
                        0,
                        0,
                        0,
                    ]
                )

                # Spawn a gamma ray emission with energy from gamma-ray list
                # energy transformed to rest frame
                cmf_energy = sample_energy(energy, intensity)

                if not cmf_energy:
                    print("No energy selected for this gamma ray!")
                    continue

                packet.nu_cmf = cmf_energy / H_CGS_KEV

                packet.nu_rf = packet.nu_cmf / doppler_gamma(
                    packet.get_direction_vector(),
                    packet.location_r,
                )

                packets.append(packet)

    return packets, energy_df_rows, energy_plot_df_rows


def main_gamma_ray_loop(num_decays, model, plasma):
    """Main loop that determines the gamma ray propagation

    Parameters
    ----------
    num_decays : int
        Number of decays requested
    model : tardis.Radial1DModel
        The tardis model to calculate gamma ray propagation through
    plasma : tardis.plasma.BasePlasma
        The tardis plasma with calculated atomic number density

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

    new_abundance_rows = []

    electron_number_density = (
        plasma.number_density.mul(plasma.number_density.index, axis=0)
    ).sum()

    # change to .decay() of raw_isotope_abundance
    for index, mass in enumerate(shell_masses):
        isotope_abundance = raw_isotope_abundance[index]
        isotope_dict = {}
        for (
            atom_number,
            atom_mass,
        ), abundance in isotope_abundance.iteritems():
            isotope_string = atomic_number2element_symbol(atom_number) + str(
                atom_mass
            )
            isotope_dict[isotope_string] = abundance
        new_abundance_rows.append(
            decay_nuclides(
                mass * u.g.to("M_sun"),
                isotope_dict,
                np.array([time_explosion * u.s.to("day")]),
            )
        )

    new_abundances = pd.concat(new_abundance_rows).transpose()
    new_abundances.columns = raw_isotope_abundance.columns

    # scale abundances, needs to be done per-isotope in future
    new_abundances *= raw_isotope_abundance.values

    (
        decays_per_shell,
        scaled_activity_df,
        activity_df,
    ) = compute_required_photons_per_shell_artis(
        shell_masses, new_abundances, num_decays
    )

    ni56_lines = read_artis_lines("ni56")
    co56_lines = read_artis_lines("co56")

    total_energy = (
        activity_df["Ni56"]
        * (ni56_lines.energy * 1000 * ni56_lines.intensity).sum()
        + activity_df["Co56"]
        * (co56_lines.energy * 1000 * co56_lines.intensity).sum()
    )

    print("Total expected energy")
    print(total_energy.sum() * u.keV.to("erg"))

    # Taking iron group to be elements 21-30
    # Used as part of the approximations for photoabsorption and pair creation
    # Dependent on atomic data
    iron_group_fraction_per_shell = model.abundance.loc[(21):(30)].sum(axis=0)

    (packets, energy_df_rows, energy_plot_df_rows,) = initialize_packets(
        number_of_shells,
        decays_per_shell,
        ejecta_volume,
        inner_velocities,
        outer_velocities,
        ni56_lines,
        co56_lines,
        scaled_activity_df,
    )

    total_cmf_energy = 0
    total_rf_energy = 0

    for p in packets:
        total_cmf_energy += p.energy_cmf
        total_rf_energy += p.energy_rf

    energy_ratio = total_energy.sum() / total_cmf_energy

    print("Total CMF energy")
    print(total_cmf_energy)

    print("Energy ratio")
    print(energy_ratio)

    for p in packets:
        p.energy_cmf *= energy_ratio
        p.energy_rf *= energy_ratio

    for e in energy_df_rows:
        e *= energy_ratio

    for row in energy_plot_df_rows:
        row[1] *= energy_ratio

    print("Total RF energy")
    print(total_rf_energy)

    i = 0
    for packet in tqdm(packets):

        while packet.status == GXPacketStatus.IN_PROCESS:

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

            # convert opacities to rest frame
            total_opacity = (
                compton_opacity
                + photoabsorption_opacity
                + pair_creation_opacity
            )

            (distance_interaction, distance_boundary,) = distance_trace(
                packet,
                inner_velocities,
                outer_velocities,
                total_opacity,
                time_explosion,
            )
            if distance_interaction < distance_boundary:

                packet.tau = -np.log(np.random.random())

                packet.status = scatter_type(
                    compton_opacity,
                    photoabsorption_opacity,
                    total_opacity,
                )

                packet = move_packet(packet, distance_interaction)
                packet.shell = np.searchsorted(
                    outer_velocities, packet.location_r, side="left"
                )
                packet.time_current += (
                    distance_interaction / C_CGS * time_explosion
                )

                packet, ejecta_energy_gained = process_packet_path(packet)

                # Save packets to dataframe rows
                # convert KeV to eV / s / cm^3
                energy_df_rows[packet.shell] += (
                    ejecta_energy_gained * 1000 / ejecta_volume[packet.shell]
                )
                energy_plot_df_rows.append(
                    [
                        i,
                        ejecta_energy_gained
                        * 1000
                        / ejecta_volume[packet.shell],
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
                packet.tau -= total_opacity * distance_boundary * time_explosion
                # overshoot so that the gamma-ray is comfortably in the next shell
                packet = move_packet(
                    packet, distance_boundary * (1 + BOUNDARY_THRESHOLD)
                )
                packet.time_current += (
                    distance_boundary / C_CGS * time_explosion
                )
                packet.shell = np.searchsorted(
                    outer_velocities, packet.location_r, side="left"
                )

            if packet.shell > len(ejecta_density) - 1:
                escape_energy.append(packet.energy_rf)
                packet.status = GXPacketStatus.END
            elif packet.shell < 0:
                packet.energy_rf = 0.0
                packet.energy_cmf = 0.0
                packet.status = GXPacketStatus.END

        i += 1

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

    # Energy is eV/s/cm^-3
    energy_df = pd.DataFrame(data=energy_df_rows, columns=["energy"])

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
        scaled_activity_df,
        activity_df,
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
