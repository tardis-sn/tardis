import astropy.units as u
import tardis.constants as const
import numpy as np
from astropy.coordinates import spherical_to_cartesian

R_ELECTRON = const.a0.cgs * const.alpha.cgs ** 2.0
KEV2ERG = (1000 * u.eV).to("erg").value


class SphericalVector(object):
    """
    Direction object to hold spherical polar and Cartesian directions
    Must be initialized with r, theta, phi

    Attributes
    ----------
    r : float64
        vector r position
    theta : float64
        vector mu position
    phi : float64
        vector phi position
    mu : float64
        calculated vector theta position
    x : float64
        calculated vector x position
    y : float64
        calculated vector y position
    z : float64
        calculated vector z position
    """

    def __init__(self, r, theta, phi=0.0):
        self.r = r
        self.theta = theta
        self.phi = phi

    @property
    def mu(self):
        return np.cos(self.theta)

    @property
    def cartesian_coords(self):
        # 0.5*np.pi subtracted because of the definition of theta
        # in astropy.coordinates.cartesian_to_spherical
        x, y, z = spherical_to_cartesian(
            self.r, self.theta - 0.5 * np.pi, self.phi
        )
        return x.value, y.value, z.value


def kappa_calculation(energy):
    """
    Calculates kappa for various other calculations
    i.e. energy normalized to electron rest energy
    511.0 KeV

    Parameters
    ----------
    energy : float

    Returns
    -------
    kappa : float

    """
    k = energy / 511.0
    return k


def euler_rodrigues(theta, direction):
    """
    Calculates the Euler-Rodrigues rotation matrix

    Parameters
    ----------
    theta : float
    direction : SphericalVector object

    Returns
    -------
    rotation matrix : Two dimensional Numpy array, dtype float

    """
    a = np.cos(theta / 2)
    dir_x, dir_y, dir_z = direction
    b = dir_x * np.sin(theta / 2)
    c = dir_y * np.sin(theta / 2)
    d = dir_z * np.sin(theta / 2)

    er11 = a ** 2.0 + b ** 2.0 - c ** 2.0 - d ** 2.0
    er12 = 2.0 * (b * c - a * d)
    er13 = 2.0 * (b * d + a * c)

    er21 = 2.0 * (b * c + a * d)
    er22 = a ** 2.0 + c ** 2.0 - b ** 2.0 - d ** 2.0
    er23 = 2.0 * (c * d - a * b)

    er31 = 2.0 * (b * d - a * c)
    er32 = 2.0 * (c * d + a * b)
    er33 = a ** 2.0 + d ** 2.0 - b ** 2.0 - c ** 2.0

    return np.array(
        [[er11, er12, er13], [er21, er22, er23], [er31, er32, er33]]
    )


def solve_quadratic_equation(x, y, z, x_dir, y_dir, z_dir, radius_velocity):
    """
    Solves the quadratic equation for the distance to the shell boundary

    Parameters
    ----------
    x,y,z : float
    x_dir, y_dir, z_dir : float
    radius_velocity : float

    Returns
    -------
    solution_1 : float
    solution_2 : float

    """
    b = 2.0 * (x * x_dir + y * y_dir + z * z_dir)
    c = -(radius_velocity ** 2) + x ** 2 + y ** 2 + z ** 2
    root = b ** 2 - 4 * c
    solution_1 = -np.inf
    solution_2 = -np.inf
    if root > 0.0:
        solution_1 = 0.5 * (-b + np.sqrt(root))
        solution_2 = 0.5 * (-b - np.sqrt(root))
    elif root == 0:
        solution_1 = -0.5 * b
    return solution_1, solution_2


def klein_nishina(energy, theta_C):
    """
    Calculates the Klein-Nishina equation

    Parameters
    ----------
    energy : float
        Packet energy
    theta_C : float
        Compton angle
    Returns
    -------
    float
        Klein-Nishina solution
    """
    kappa = kappa_calculation(energy)
    return (
        R_ELECTRON.value
        / 2
        * (1.0 + kappa * (1.0 - np.cos(theta_C))) ** -2.0
        * (
            1.0
            + np.cos(theta_C) ** 2.0
            + (kappa ** 2.0 * (1.0 - np.cos(theta_C)) ** 2.0)
            / (1.0 + kappa * (1.0 - np.cos(theta_C)))
        )
    )


def compton_theta_distribution(energy, sample_resolution=100):
    """
    Calculates the cumulative distribution function of theta angles
    for Compton Scattering

    Parameters
    ----------
    energy : float
    sample_resolution : int

    Returns
    -------
    theta_angles : One dimensional Numpy array, dtype float
    norm_theta_distribution : One dimensional Numpy array, dtype float

    """
    dtheta = np.pi / sample_resolution

    theta_distribution = np.zeros(sample_resolution)
    theta_angles = np.ones(sample_resolution) * np.pi

    theta = 0.0
    for i in range(sample_resolution - 1):
        theta_distribution[i + 1] = theta_distribution[i] + klein_nishina(
            energy, theta
        )
        theta_angles[i] = theta
        theta += dtheta

    norm_theta_distribution = theta_distribution / np.max(theta_distribution)

    return theta_angles, norm_theta_distribution


def get_random_theta_gamma_ray():
    """Get a random theta direction between 0 and pi
    Returns
    -------
    float
        Random theta direction
    """
    return np.arccos(1.0 - 2.0 * np.random.random())


def get_random_phi_gamma_ray():
    """Get a random phi direction between 0 and 2 * pi

    Returns
    -------
    float
        Random phi direction
    """
    return 2.0 * np.pi * np.random.random()


def convert_half_life_to_astropy_units(half_life_string):
    """Converts input half-life to use astropy units

    Parameters
    ----------
    half_life_string : string
        Half-life as a string

    Returns
    -------
    astropy unit
        Half-life in seconds
    """
    [value, unit] = half_life_string.split(" ")
    try:
        [nominal_value, magnitude] = value.split("×")
        [base, exponent] = magnitude.split("+")
        nominal_value = float(nominal_value) * float(base) ** float(exponent)
    except ValueError:
        nominal_value = float(value)
    if unit == "y":
        unit = "yr"
    half_life_with_unit = nominal_value * u.Unit(unit)
    return half_life_with_unit.to(u.s)


def calculate_energy_per_mass(energy, mass):
    """Calculates energy in erg / g

    Parameters
    ----------
    energy : float64
        Energy in KeV
    mass : float64
        Mass in grams

    Returns
    -------
    float64
        Energy density in erg / g
    """
    return energy * KEV2ERG / mass


def normalize(vector):
    """
    Normalizes a vector in cartesian coordinates

    Parameters
    ----------
    vector : One-dimensional Numpy Array, dtype float

    Returns
    -------
    normalized_vector : One-dimensional Numpy Array, dtype float
    """
    normalized_vector = vector / np.sqrt(
        vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2
    )
    return normalized_vector


def get_perpendicular_vector(original_direction):
    """
    Computes a vector which is perpendicular to the input vector

    Parameters
    ----------
    original_direction : SphericalVector object

    Returns
    -------

    """
    # draw random angles
    theta = get_random_theta_gamma_ray()
    phi = get_random_phi_gamma_ray()
    # transform random angles to cartesian coordinates
    # 0.5*np.pi subtracted because of the definition of theta
    # in astropy.coordinates.cartesian_to_spherical
    random_vector = spherical_to_cartesian(1, theta - 0.5 * np.pi, phi)
    perpendicular_vector = np.cross(original_direction, random_vector)
    perpendicular_vector = normalize(perpendicular_vector)
    return perpendicular_vector
