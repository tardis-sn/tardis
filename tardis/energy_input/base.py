import numpy as np

from tardis.montecarlo.montecarlo_numba.r_packet import get_random_mu
from tardis.energy_input.util import SphericalVector
from tardis.energy_input.gamma_ray_grid import density_sampler, distance_trace, move_gamma_ray, mass_distribution
from tardis.energy_input.energy_source import setup_gamma_ray_energy, read_nuclear_dataframe, sample_energy_distribution
from tardis.energy_input.calculate_opacity import compton_opacity_calculation, photoabsorption_opacity_calculation, pair_creation_opacity_calculation, kappa_calculation
from tardis.energy_input.gamma_ray_interactions import scatter_type
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

def spawn_gamma_ray(gamma_ray, radii, mass_ratio, energy_sorted, energy_cdf):
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
    gamma_ray : GammaRay object

    """

    direction_mu = get_random_mu()
    direction_phi = 0.0
    
    energy_KeV = sample_energy_distribution(energy_sorted, energy_cdf)

    gamma_ray.energy = energy_KeV
        
    initial_radius, shell = density_sampler(radii, mass_ratio)
    
    if shell < len(radii) - 1:
        initial_radius += np.random.random() * (radii[shell + 1] - radii[shell])
    
    gamma_ray.shell = shell
        
    gamma_ray.direction = SphericalVector(1., direction_mu, direction_phi)
    
    location_mu = get_random_mu()
    location_phi = 0.0
    gamma_ray.location = SphericalVector(initial_radius, location_mu, location_phi)
        
    return gamma_ray


def main_gamma_ray_loop(num_packets, model, path):

    output_energies = []
    ejecta_energy = []
    ejecta_energy_r = []

    inner_radius = model.r_inner[0].value
    outer_radius =  model.r_outer[-1].value
    
    outer_radii = model.r_outer[:].value
    inner_radii = model.r_inner[:].value

    ejecta_density = model.density[:].value

    masses = mass_distribution(model.no_of_shells, inner_radii, ejecta_density)

    iron_group_fraction = 0.5

    packets = []

    nuclear_data = read_nuclear_dataframe(path)

    energy_sorted, energy_cdf = setup_gamma_ray_energy(nuclear_data)

    for i in range(num_packets):
        
        ray = GammaRay(0, 0, 1, 'InProcess', 0)
        packets.append(spawn_gamma_ray(ray, outer_radii, masses, energy_sorted, energy_cdf))

    i=0
    for packet in tqdm(packets):

        distance_moved = 0.

        while packet.status == "InProcess":
            compton_opacity = compton_opacity_calculation(ejecta_density[packet.shell], packet.energy)
            photoabsorption_opacity = photoabsorption_opacity_calculation(packet.energy, ejecta_density[packet.shell], iron_group_fraction)
            pair_creation_opacity = pair_creation_opacity_calculation(packet.energy, ejecta_density[packet.shell], iron_group_fraction)
            total_opacity = compton_opacity + photoabsorption_opacity + pair_creation_opacity
            
            distance_interaction, distance_boundary, interaction = \
                                            distance_trace(packet, inner_radii, outer_radii, total_opacity, distance_moved)
            
            if interaction:
                ejecta_energy_gained, pair_created = scatter_type(packet, compton_opacity, photoabsorption_opacity, total_opacity)
                #Add antiparallel packet on pair creation at end of list
                if pair_created:
                    backward_ray = packet
                    backward_ray.direction.phi += np.pi
                    packets.append(backward_ray)

                if ejecta_energy_gained > 0.0:
                    ejecta_energy.append(ejecta_energy_gained)
                    ejecta_energy_r.append(packet.location.r)
                    
                packet = move_gamma_ray(packet, distance_interaction)
                distance_moved = 0.
                
            else:
                rad_before = packet.location.r
                packet = move_gamma_ray(packet, distance_boundary)
                rad_after = packet.location.r
                distance_moved = distance_boundary
                if rad_after > rad_before:
                    packet.shell += 1
                else:
                    packet.shell -= 1
                    
            if (packet.location.r - outer_radius) < 1.0 or packet.shell > len(ejecta_density) - 1:
                packet.status = 'Emitted'
                output_energies.append(packet.energy)
            elif packet.location.r < inner_radius:
                packet.status = 'Absorbed'
                packet.energy = 0.0
            
            if packet.status == 'PhotoAbsorbed':
                #log where energy is deposited
                ejecta_energy.append(packet.energy)
                ejecta_energy_r.append(packet.location.r)
        
        i+=1

    return ejecta_energy, ejecta_energy_r, output_energies, inner_radii