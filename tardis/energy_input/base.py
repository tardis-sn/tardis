import numpy as np
import copy
from tqdm.auto import tqdm
import pandas as pd

# from tardis.montecarlo.montecarlo_numba.r_packet import get_random_mu
from tardis.energy_input.util import SphericalVector
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
)
from tardis.energy_input.util import (
    get_random_theta_gamma_ray,
    get_random_phi_gamma_ray,
    calculate_energy_per_mass,
)
from tardis import constants as const
from astropy.coordinates import cartesian_to_spherical


class GXPacket(object):
    """
    Gamma ray or X ray object with location, direction, energy, time and optical depth

    Attributes
    ----------
    location : SphericalVector object
             GXPacket position vector
    direction : SphericalVector object
             GXPacket direction vector (unitary)
    energy : float64
             GXPacket energy
    status : str
             GXPacket status
    shell : int64
             GXPacket shell location index
    """

    def __init__(self, location, direction, energy, status, shell):
        self.location = location
        self.direction = direction
        self.energy = energy
        self.status = status
        self.shell = shell
        self.time_created = 0
        self.time_current = 0
        self.tau = -np.log(np.random.random())


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
        Energy input per packet,
        radius of deposition,
        theta angle of deposition,
        time of deposition,
        type of deposition where:
            -1 = beta decay,
            0 = Compton scatter,
            1 = photoabsorption,
            2 = pair creation

    """
    escape_energy = []
    interaction_count = []

    # Note the use of velocity as the radial coordinate
    inner_radius = model.v_inner[0].value
    outer_radius = model.v_outer[-1].value
    outer_radii = model.v_outer[:].value
    inner_radii = model.v_inner[:].value
    ejecta_density = model.density[:].value
    ejecta_epoch = model.time_explosion.to("s").value
    number_of_shells = model.no_of_shells
    raw_isotope_abundance = model.raw_isotope_abundance

    shell_masses = mass_per_shell(
        number_of_shells, inner_radii, outer_radii, ejecta_density
    )

    packets_per_shell, decay_rad_db = compute_required_packets_per_shell(
        outer_radii,
        inner_radii,
        ejecta_density,
        number_of_shells,
        raw_isotope_abundance,
        num_packets,
    )

    # Taking iron group to be elements 21-30
    iron_group_fraction_per_shell = raw_isotope_abundance.loc[(21,):(30,)].sum(
        axis=0
    )

    packets = []
    energy_plot_df_rows = []
    energy_df_rows = np.zeros(number_of_shells)

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
        for shell in range(model.no_of_shells):
            required_packets_in_shell = packets_per_shell[column].iloc[shell]
            for i_packet in range(required_packets_in_shell):
                # draw a random gamma-ray in shell
                ray1 = GXPacket(0, 0, 1, "InProcess", 0)
                z1 = np.random.random()
                initial_radius = inner_radii[shell] + z1 * (
                    outer_radii[shell] - inner_radii[shell]
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
                    ) = ray1.direction.get_cartesian_coords
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

    i = 0
    for packet in tqdm(packets):
        distance_moved = 0.0
        interaction_count.append(0)

        while packet.status == "InProcess":

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
                inner_radii,
                outer_radii,
                total_opacity,
                ejecta_epoch,
            )

            if distance_interaction < distance_boundary:
                interaction_count[i] += 1

                packet.tau = -np.log(np.random.random())

                ejecta_energy_gained, compton_angle = scatter_type(
                    packet,
                    compton_opacity,
                    photoabsorption_opacity,
                    total_opacity,
                )

                if (
                    packet.status == "ComptonScatter"
                    and ejecta_energy_gained > 0.0
                ):
                    packet = move_gamma_ray(packet, distance_interaction)
                    packet.shell = get_shell(packet.location.r, outer_radii)
                    distance_moved += distance_interaction * ejecta_epoch
                    packet.time_current += (
                        distance_interaction / const.c.cgs.value * ejecta_epoch
                    )

                    (
                        packet.direction.theta,
                        packet.direction.phi,
                    ) = compton_scatter(packet, compton_angle)

                    energy_df_rows[packet.shell] += calculate_energy_per_mass(
                        ejecta_energy_gained, shell_masses[shell]
                    )
                    energy_plot_df_rows.append(
                        [
                            i,
                            ejecta_energy_gained,
                            packet.location.r,
                            packet.location.theta,
                            packet.time_current,
                            0,
                        ]
                    )

                if (
                    packet.status == "PhotoAbsorbed"
                    and ejecta_energy_gained > 0.0
                ):
                    packet = move_gamma_ray(packet, distance_interaction)
                    packet.shell = get_shell(packet.location.r, outer_radii)
                    distance_moved += distance_interaction * ejecta_epoch
                    packet.time_current += (
                        distance_interaction / const.c.cgs.value * ejecta_epoch
                    )

                    energy_df_rows[packet.shell] += calculate_energy_per_mass(
                        ejecta_energy_gained, shell_masses[shell]
                    )
                    energy_plot_df_rows.append(
                        [
                            i,
                            ejecta_energy_gained,
                            packet.location.r,
                            packet.location.theta,
                            packet.time_current,
                            1,
                        ]
                    )

                    # Packet destroyed, go to the next packet

                if packet.status == "PairCreated":
                    packet = move_gamma_ray(packet, distance_interaction)
                    packet.shell = get_shell(packet.location.r, outer_radii)
                    distance_moved += distance_interaction * ejecta_epoch
                    packet.time_current += (
                        distance_interaction / const.c.cgs.value * ejecta_epoch
                    )

                    energy_df_rows[packet.shell] += calculate_energy_per_mass(
                        ejecta_energy_gained, shell_masses[shell]
                    )
                    energy_plot_df_rows.append(
                        [
                            i,
                            ejecta_energy_gained,
                            packet.location.r,
                            packet.location.theta,
                            packet.time_current,
                            2,
                        ]
                    )

                    pair_creation(packet)

                    backward_ray = GXPacket(
                        copy.deepcopy(packet.location),
                        copy.deepcopy(packet.direction),
                        copy.deepcopy(packet.energy),
                        "InProcess",
                        copy.deepcopy(packet.shell),
                    )

                    backward_ray.direction.phi += np.pi

                    if backward_ray.direction.phi > 2 * np.pi:
                        backward_ray.direction.phi -= 2 * np.pi

                    # Add antiparallel packet on pair creation at end of list
                    packets.append(backward_ray)

                if packet.status == "PhotoAbsorbed":
                    break
                else:
                    packet.status = "InProcess"

                distance_moved = 0.0

            else:
                packet.tau -= total_opacity * distance_boundary * ejecta_epoch
                packet = move_gamma_ray(packet, distance_boundary)
                distance_moved += distance_boundary * ejecta_epoch
                packet.time_current += (
                    distance_boundary / const.c.cgs.value * ejecta_epoch
                )
                packet.shell = get_shell(packet.location.r, outer_radii)
                if distance_boundary == 0.0:
                    packet.shell += 1

            if (
                np.abs(packet.location.r - outer_radius) < 10.0
                or packet.shell > len(ejecta_density) - 1
            ):
                packet.status = "Emitted"
                escape_energy.append(packet.energy)
            elif (
                np.abs(packet.location.r - inner_radius) < 10.0
                or packet.shell < 0
            ):
                packet.status = "Absorbed"
                packet.energy = 0.0

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

    return (energy_df, energy_plot_df)
