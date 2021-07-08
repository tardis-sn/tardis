import numpy as np
import copy
from tqdm.auto import tqdm
import pandas as pd

# from tardis.montecarlo.montecarlo_numba.r_packet import get_random_mu
from tardis.energy_input.util import SphericalVector, GXPacket, GXPacketStatus
from tardis.energy_input.gamma_ray_grid import (
    distance_trace,
    move_gamma_ray,
    get_shell,
    compute_required_packets_per_shell,
    mass_per_shell,
)
from tardis.energy_input.energy_source import (
    setup_input_energy,
    sample_energy_distribution,
    intensity_ratio,
)
from tardis.energy_input.calculate_opacity import (
    compton_opacity_calculation,
    photoabsorption_opacity_calculation,
    pair_creation_opacity_calculation,
)
from tardis.energy_input.gamma_ray_interactions import (
    scatter_type,
    compton_scatter,
    pair_creation,
    get_compton_angle,
)
from tardis.energy_input.util import (
    get_random_theta_gamma_ray,
    get_random_phi_gamma_ray,
    calculate_energy_per_mass,
)
from tardis import constants as const
from astropy.coordinates import cartesian_to_spherical
from tardis.montecarlo.montecarlo_numba.numba_config import CLOSE_LINE_THRESHOLD


def initialize_packets(
    number_of_shells,
    packets_per_shell,
    shell_masses,
    inner_velocities,
    outer_velocities,
    decay_rad_db,
):
    """Initializes packet properties
    and appends beta decay energy to output tables

    Parameters
    ----------
    number_of_shells : int64
        Number of shells in model
    packets_per_shell : pandas dataframe
        Number of packets in a shell
    shell_masses : array of float64
        Mass per shell
    inner_velocities : array of float64
        Shell inner velocities
    outer_velocities : array of float64
        Shell outer velocities
    decay_rad_db : pandas dataframe
        Decay radiation database

    Returns
    -------
    list
        GXPacket objects
    numpy array
        energy binned per shell
    list
        packet info
    """
    packets = []
    energy_df_rows = np.zeros(number_of_shells)
    energy_plot_df_rows = []

    for column in packets_per_shell:
        subtable = decay_rad_db.loc[column]
        gamma_ratio, positron_ratio = intensity_ratio(
            subtable, "'gamma_rays' or type=='x_rays'", "'e+'"
        )
        energy_sorted, energy_cdf = setup_input_energy(
            subtable, "'gamma_rays' or type=='x_rays'"
        )
        positron_energy_sorted, positron_energy_cdf = setup_input_energy(
            subtable, "'e+'"
        )
        for shell in range(number_of_shells):
            required_packets_in_shell = packets_per_shell[column].iloc[shell]
            for i_packet in range(required_packets_in_shell):
                # draw a random gamma-ray in shell
                ray1 = GXPacket(0, 0, 1, GXPacketStatus.IN_PROCESS, 0)
                z1 = np.random.random()
                initial_radius = inner_velocities[shell] + z1 * (
                    outer_velocities[shell] - inner_velocities[shell]
                )
                location_theta = get_random_theta_gamma_ray()
                location_phi = get_random_phi_gamma_ray()
                ray1.location = SphericalVector(
                    initial_radius, location_theta, location_phi
                )
                direction_theta = get_random_theta_gamma_ray()
                direction_phi = get_random_phi_gamma_ray()
                ray1.direction = SphericalVector(
                    1.0, direction_theta, direction_phi
                )
                ray1.shell = shell
                z2 = np.random.random()
                if z2 > gamma_ratio:
                    # positron: sets gamma-ray energy to 511keV
                    ray1.energy = 511.0
                    packets.append(ray1)

                    # annihilation dumps energy into medium
                    energy_KeV = sample_energy_distribution(
                        positron_energy_sorted, positron_energy_cdf
                    )

                    energy_df_rows[shell] += calculate_energy_per_mass(
                        energy_KeV, shell_masses[shell]
                    )
                    energy_plot_df_rows.append(
                        [
                            -1,
                            energy_KeV,
                            initial_radius,
                            ray1.location.theta,
                            0.0,
                            -1,
                        ]
                    )

                    # annihilation produces second gamma-ray in opposite direction
                    (
                        ray1_x,
                        ray1_y,
                        ray1_z,
                    ) = ray1.direction.cartesian_coords
                    ray2 = copy.deepcopy(ray1)
                    ray2_x = -ray1_x
                    ray2_y = -ray1_y
                    ray2_z = -ray1_z
                    ray2_r, ray2_theta, ray2_phi = cartesian_to_spherical(
                        ray2_x, ray2_y, ray2_z
                    )
                    ray2.direction.r = ray2_r
                    ray2.direction.theta = ray2_theta.value + 0.5 * np.pi
                    ray2.direction.phi = ray2_phi.value

                else:
                    # Spawn a gamma ray emission with energy from gamma-ray list
                    ray1.energy = sample_energy_distribution(
                        energy_sorted, energy_cdf
                    )
                    packets.append(ray1)

    return packets, energy_df_rows, energy_plot_df_rows


def main_gamma_ray_loop(num_packets, model):
    """Main loop that determines the gamma ray propagation

    Parameters
    ----------
    num_packets : int
        Number of packets
    model : tardis.Radial1DModel
        The tardis model to calculate gamma ray propagation through

    Returns
    -------
    pandas DataFrame
        Energy per mass per shell in units of erg / g
    pandas DataFrame
        Columns:
        Packet index,
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
    """
    escape_energy = []
    interaction_count = []

    # Note the use of velocity as the radial coordinate
    inner_radius = model.v_inner[0].value
    outer_radius = model.v_outer[-1].value
    outer_velocities = model.v_outer[:].value
    inner_velocities = model.v_inner[:].value
    ejecta_density = model.density[:].value
    ejecta_volume = model.volume[:].value
    ejecta_epoch = model.time_explosion.to("s").value
    number_of_shells = model.no_of_shells
    raw_isotope_abundance = model.raw_isotope_abundance

    shell_masses = mass_per_shell(ejecta_volume, ejecta_density)

    packets_per_shell, decay_rad_db = compute_required_packets_per_shell(
        shell_masses,
        raw_isotope_abundance,
        num_packets,
    )

    # Taking iron group to be elements 21-30
    iron_group_fraction_per_shell = raw_isotope_abundance.loc[(21,):(30,)].sum(
        axis=0
    )

    packets, energy_df_rows, energy_plot_df_rows = initialize_packets(
        number_of_shells,
        packets_per_shell,
        shell_masses,
        inner_velocities,
        outer_velocities,
        decay_rad_db,
    )

    i = 0
    for packet in tqdm(packets):
        interaction_count.append(0)

        while packet.status == GXPacketStatus.IN_PROCESS:

            compton_opacity = compton_opacity_calculation(
                packet.energy, ejecta_density[packet.shell]
            )
            photoabsorption_opacity = photoabsorption_opacity_calculation(
                packet.energy,
                ejecta_density[packet.shell],
                iron_group_fraction_per_shell[packet.shell],
            )
            pair_creation_opacity = pair_creation_opacity_calculation(
                packet.energy,
                ejecta_density[packet.shell],
                iron_group_fraction_per_shell[packet.shell],
            )
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
                ejecta_epoch,
            )

            if distance_interaction < distance_boundary:
                interaction_count[i] += 1

                packet.tau = -np.log(np.random.random())

                packet.status = scatter_type(
                    compton_opacity,
                    photoabsorption_opacity,
                    total_opacity,
                )

                packet = move_gamma_ray(packet, distance_interaction)
                packet.shell = get_shell(packet.location.r, outer_velocities)
                packet.time_current += (
                    distance_interaction / const.c.cgs.value * ejecta_epoch
                )

                if packet.status == GXPacketStatus.COMPTON_SCATTER:
                    (
                        compton_angle,
                        ejecta_energy_gained,
                        packet.energy,
                    ) = get_compton_angle(packet.energy)
                    (
                        packet.direction.theta,
                        packet.direction.phi,
                    ) = compton_scatter(packet, compton_angle)

                if packet.status == GXPacketStatus.PAIR_CREATION:
                    ejecta_energy_gained = packet.energy - (2.0 * 511.0)
                    packet, backward_ray = pair_creation(packet)

                    # Add antiparallel packet on pair creation at end of list
                    packets.append(backward_ray)

                if packet.status == GXPacketStatus.PHOTOABSORPTION:
                    ejecta_energy_gained = packet.energy

                # Save packets to dataframe rows
                energy_df_rows[packet.shell] += calculate_energy_per_mass(
                    ejecta_energy_gained, shell_masses[packet.shell]
                )
                energy_plot_df_rows.append(
                    [
                        i,
                        ejecta_energy_gained,
                        packet.location.r,
                        packet.location.theta,
                        packet.time_current,
                        int(packet.status),
                    ]
                )

                if packet.status == GXPacketStatus.PHOTOABSORPTION:
                    # Packet destroyed, go to the next packet
                    break
                else:
                    packet.status = GXPacketStatus.IN_PROCESS

            else:
                packet.tau -= total_opacity * distance_boundary * ejecta_epoch
                # overshoot so that the gamma-ray is comfortably in the next shell
                packet = move_gamma_ray(packet, distance_boundary * (1 + 1e-7))
                packet.time_current += (
                    distance_boundary / const.c.cgs.value * ejecta_epoch
                )
                packet.shell = get_shell(packet.location.r, outer_velocities)

            if (
                np.abs(packet.location.r - outer_radius) < 10.0
                or packet.shell > len(ejecta_density) - 1
            ):
                escape_energy.append(packet.energy)
                packet.status = GXPacketStatus.END
            elif (
                np.abs(packet.location.r - inner_radius) < 10.0
                or packet.shell < 0
            ):
                packet.energy = 0.0
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
        ],
    )

    energy_df = pd.DataFrame(data=energy_df_rows, columns=["energy_per_mass"])

    return (energy_df, energy_plot_df, escape_energy)
