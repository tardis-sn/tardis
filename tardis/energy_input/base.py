import numpy as np
import copy

# from tardis.montecarlo.montecarlo_numba.r_packet import get_random_mu
from tardis.energy_input.util import SphericalVector
from tardis.energy_input.gamma_ray_grid import (
    density_sampler,
    distance_trace,
    move_gamma_ray,
    mass_distribution,
    get_shell,
)
from tardis.energy_input.energy_source import (
    setup_input_energy,
    read_nuclear_dataframe,
    sample_energy_distribution,
    intensity_ratio,
)
from tardis.energy_input.calculate_opacity import (
    compton_opacity_calculation,
    photoabsorption_opacity_calculation,
    pair_creation_opacity_calculation,
    kappa_calculation,
)
from tardis.energy_input.gamma_ray_interactions import scatter_type
from tardis.energy_input.util import get_random_mu_gamma_ray
from tqdm.auto import tqdm


class GammaRay(object):
    """
    Gamma ray object with location, direction and energy

    Attributes
    ----------
    location : SphericalVector object
             GammaRay position vector
    direction : SphericalVector object
             GammaRay direction vector (unitary)
    energy : float64
             GammaRay energy
    status : str
             GammaRay status
    shell : int64
             GammaRay shell location index
    """

    def __init__(self, location, direction, energy, status, shell):
        self.location = location
        self.direction = direction
        self.energy = energy
        self.status = status
        self.shell = shell


def spawn_gamma_ray(
    gamma_ray,
    inner_radii,
    outer_radii,
    mass_ratio,
    energy_sorted,
    energy_cdf,
    beta_decay=False,
):
    """
    Initializes a gamma ray in the simulation grid

    Parameters
    ----------
    gamma_ray : GammaRay object
    radii : One-dimensional Numpy Array, dtype float
    mass_ratio : One-dimensional Numpy Array, dtype int
    positron : dtype bool

    Returns
    -------
    GammaRay object

    """

    direction_mu = get_random_mu_gamma_ray()
    direction_phi = 0.0

    if beta_decay:
        location_mu = get_random_mu_gamma_ray()
        location_phi = 0.0
        gamma_ray.location = SphericalVector(0, location_mu, location_phi)
        gamma_ray.direction = SphericalVector(1.0, direction_mu, direction_phi)
        gamma_ray.energy = 511.0
        return gamma_ray

    energy_KeV = sample_energy_distribution(energy_sorted, energy_cdf)

    gamma_ray.energy = energy_KeV

    initial_radius, shell = density_sampler(inner_radii, mass_ratio)

    initial_radius += np.random.random() * (
        outer_radii[shell] - inner_radii[shell]
    )

    gamma_ray.shell = shell

    gamma_ray.direction = SphericalVector(1.0, direction_mu, direction_phi)

    location_mu = get_random_mu_gamma_ray()
    location_phi = 0.0
    gamma_ray.location = SphericalVector(
        initial_radius, location_mu, location_phi
    )

    return gamma_ray


def spawn_positron(
    inner_radii, outer_radii, mass_ratio, energy_sorted, energy_cdf
):
    """Spawns a positron and immediately decays it to an energy

    Parameters
    ----------
    inner_radii : float
        Inner radius of simulation
    outer_radii : float
        Outer radius of simulation
    mass_ratio : One-dimensional Numpy Array, dtype float
        CDF of masses
    energy_sorted : One-dimensional Numpy Array, dtype float
        Sorted array of possible energies
    energy_cdf : One-dimensional Numpy Array, dtype float
        CDF of possible energies

    Returns
    -------
    float
        Emitted energy
    float
        Radius of emitted energy
    int
        Shell of emitted energy
    """
    energy_KeV = sample_energy_distribution(energy_sorted, energy_cdf)

    initial_radius, shell = density_sampler(inner_radii, mass_ratio)

    initial_radius += np.random.random() * (
        outer_radii[shell] - inner_radii[shell]
    )

    return energy_KeV, initial_radius, shell


def main_gamma_ray_loop(num_packets, model, path):
    """Main loop that determines the gamma ray propagation

    Parameters
    ----------
    num_packets : int
        Number of gamma-ray packets to run
    model : tardis.Radial1DModel
        Input tardis model
    path : str
        Path to nuclear dataframe

    Returns
    -------
    list
        Energy deposited in the ejecta from events
    list
        Radial location of energy deposition
    list
        Theta location of energy deposition
    list
        Escaped energy
    One-dimensional Numpy Array, dtype float
        Inner radii of ejecta grid shells
    list
        Event type that deposited energy
    list
        Number of events each packet encountered
    list
        Initial positions of packets
    """
    output_energies = []
    ejecta_energy = []
    ejecta_energy_r = []
    ejecta_energy_theta = []
    # list of energy input types as integers 0, 1, 2 for Compton scattering, photoabsorption, pair creation
    energy_input_type = []
    interaction_count = []
    initial_positions = []

    inner_radius = model.r_inner[0].value
    outer_radius = model.r_outer[-1].value

    outer_radii = model.r_outer[:].value
    inner_radii = model.r_inner[:].value

    ejecta_density = model.density[:].value

    mass_cdf = mass_distribution(
        model.no_of_shells, inner_radii, outer_radii, ejecta_density
    )

    iron_group_fraction = 0.5

    packets = []

    nuclear_data = read_nuclear_dataframe(path)

    gamma_ratio, positron_ratio = intensity_ratio(
        nuclear_data, "'gamma_rays'", "'e+'"
    )

    # Need to decay particles, not just spawn gamma rays
    energy_sorted, energy_cdf = setup_input_energy(nuclear_data, "'gamma_rays'")
    positron_energy_sorted, positron_energy_cdf = setup_input_energy(
        nuclear_data, "'e+'"
    )

    for i in range(num_packets):

        z = np.random.random()

        if z > gamma_ratio:
            # Spawn a pair annihilation
            energy_KeV, initial_radius, shell = spawn_positron(
                inner_radii,
                outer_radii,
                mass_cdf,
                positron_energy_sorted,
                positron_energy_cdf,
            )

            # Dump energy into the ejecta (via presumed Coulomb interaction)
            energy_input_type.append(2)
            ejecta_energy.append(energy_KeV)
            ejecta_energy_r.append(initial_radius)

            ray = GammaRay(0, 0, 1, "InProcess", 0)
            gamma_ray = spawn_gamma_ray(
                ray,
                inner_radii,
                outer_radii,
                mass_cdf,
                energy_sorted,
                energy_cdf,
                beta_decay=True,
            )

            gamma_ray.location.r = initial_radius
            gamma_ray.shell = shell

            ejecta_energy_theta.append(gamma_ray.location.theta)

            packets.append(gamma_ray)

            gamma_ray_2 = copy.deepcopy(gamma_ray)
            gamma_ray_2.direction.phi += np.pi

            if gamma_ray_2.direction.phi > 2 * np.pi:
                gamma_ray_2.direction.phi -= 2 * np.pi

            packets.append(gamma_ray_2)

        else:
            # Spawn a gamma ray emission
            ray = GammaRay(0, 0, 1, "InProcess", 0)
            gamma_ray = spawn_gamma_ray(
                ray,
                inner_radii,
                outer_radii,
                mass_cdf,
                energy_sorted,
                energy_cdf,
            )
            packets.append(gamma_ray)

    i = 0
    for packet in tqdm(packets):
        initial_positions.append(packet.location.r)

        distance_moved = 0.0

        z = np.random.random()
        tau = -np.log(z)

        interaction_count.append(0)
        j = 0
        while packet.status == "InProcess":

            compton_opacity = compton_opacity_calculation(
                packet.energy, ejecta_density[packet.shell]
            )
            photoabsorption_opacity = photoabsorption_opacity_calculation(
                packet.energy, ejecta_density[packet.shell], iron_group_fraction
            )
            pair_creation_opacity = pair_creation_opacity_calculation(
                packet.energy, ejecta_density[packet.shell], iron_group_fraction
            )
            total_opacity = (
                compton_opacity
                + photoabsorption_opacity
                + pair_creation_opacity
            )

            (
                distance_interaction,
                distance_boundary,
                interaction,
            ) = distance_trace(
                packet,
                tau,
                inner_radii,
                outer_radii,
                total_opacity,
                distance_moved,
            )

            # if j > 20:
            # print("Packet", i, "going slow")
            # print(np.rad2deg(packet.direction.phi))
            # print(np.rad2deg(packet.direction.theta))
            # print(distance_boundary)

            if interaction:
                interaction_count[i] += 1

                z = np.random.random()
                tau = -np.log(z)

                ejecta_energy_gained = scatter_type(
                    packet,
                    compton_opacity,
                    photoabsorption_opacity,
                    total_opacity,
                )
                # Add antiparallel packet on pair creation at end of list
                if (
                    packet.status == "ComptonScatter"
                    and ejecta_energy_gained > 0.0
                ):
                    energy_input_type.append(0)
                    packet = move_gamma_ray(packet, distance_interaction)
                    ejecta_energy.append(ejecta_energy_gained)
                    ejecta_energy_r.append(packet.location.r)
                    ejecta_energy_theta.append(packet.location.theta)

                if (
                    packet.status == "PhotoAbsorbed"
                    and ejecta_energy_gained > 0.0
                ):
                    energy_input_type.append(1)
                    packet = move_gamma_ray(packet, distance_interaction)
                    ejecta_energy.append(ejecta_energy_gained)
                    ejecta_energy_r.append(packet.location.r)
                    ejecta_energy_theta.append(packet.location.theta)
                    # Packet destroyed, go to the next packet
                    break

                if packet.status == "PairCreated":
                    packet = move_gamma_ray(packet, distance_interaction)
                    backward_ray = GammaRay(
                        packet.location,
                        packet.direction,
                        packet.energy,
                        "InProcess",
                        packet.shell,
                    )
                    backward_ray.direction.phi += np.pi
                    if backward_ray.direction.phi > 2 * np.pi:
                        backward_ray.direction.phi -= 2 * np.pi

                    packets.append(backward_ray)

                packet.status = "InProcess"
                distance_moved = 0.0

            else:
                packet = move_gamma_ray(packet, distance_boundary)
                distance_moved += distance_boundary
                packet.shell = get_shell(packet.location.r, outer_radii)

                if distance_boundary == 0.0:
                    packet.shell += 1

            if (
                np.abs(packet.location.r - outer_radius) < 10.0
                or packet.shell > len(ejecta_density) - 1
            ):
                packet.status = "Emitted"
                output_energies.append(packet.energy)
            elif (
                np.abs(packet.location.r - inner_radius) < 10.0
                or packet.shell < 0
            ):
                packet.status = "Absorbed"
                packet.energy = 0.0

            j += 1

        i += 1

    return (
        ejecta_energy,
        ejecta_energy_r,
        ejecta_energy_theta,
        output_energies,
        inner_radii,
        energy_input_type,
        interaction_count,
        initial_positions,
    )
