import numpy as np

#from tardis.montecarlo.montecarlo_numba.r_packet import get_random_mu
from tardis.energy_input.util import get_random_mu
from tardis.energy_input.gamma_ray_grid import density_sampler
from tardis.energy_input.energy_source import sample_energy_distribution

class SphericalVector(object):
    """
    Direction object to hold spherical polar and Cartesian directions
    Must be initialized with r, mu, phi

    Attributes
    ----------
    r : float64
             vector r position
    mu : float64
             vector mu position
    phi : float64
             vector phi position
    theta : float64
             calculated vector theta position
    x : float64
             calculated vector x position
    y : float64
             calculated vector y position
    z : float64
             calculated vector z position
    """
    def __init__(self, r, mu, phi=0.):
        self.r = r
        self.mu = mu
        self.phi = phi
        
    @property
    def theta(self):
        return np.arccos(self.mu)
        
    @property
    def x(self):
        return self.r * np.sin(np.arccos(self.mu)) * np.cos(self.phi)
    
    @property
    def y(self):
        return self.r * np.sin(np.arccos(self.mu)) * np.sin(self.phi)
    
    @property
    def z(self):
        return self.r * self.mu

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

def spawn_gamma_ray(gamma_ray, radii, mass_ratio, positron=False):
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
    
    if positron:
        gamma_ray.energy = 511.0e3
    else:
        gamma_ray.energy = sample_energy_distribution()
        
    initial_radius, shell = density_sampler(radii, mass_ratio)
    
    if shell < len(radii) - 1:
        initial_radius += np.random.random() * (radii[shell + 1] - radii[shell])
        
    
    gamma_ray.shell = shell
        
    gamma_ray.direction = SphericalVector(1., direction_mu, direction_phi)
    
    location_mu = get_random_mu()
    location_phi = 0.0
    gamma_ray.location = SphericalVector(initial_radius, location_mu, location_phi)
        
    return gamma_ray