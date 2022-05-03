import astropy.units as u
import tardis.constants as const
import numpy as np
from numba import njit

R_ELECTRON_SQUARED = (const.a0.cgs.value * const.alpha.cgs.value ** 2.0) ** 2.0
ELECTRON_MASS_ENERGY_KEV = (const.m_e * const.c ** 2.0).to("keV").value
BOUNDARY_THRESHOLD = 1e-7
KEV2ERG = (1000 * u.eV).to("erg").value
C_CGS = const.c.cgs.value
H_CGS_KEV = const.h.to("keV s").value


@njit
def spherical_to_cartesian(r, theta, phi):
    """Converts spherical coordinates to Cartesian

    Parameters
    ----------
    r : float64
        Radius
    theta : float64
        Theta angle in radians
    phi : float64
        Phi angle in radians

    Returns
    -------
    float64, float64, float64
        x, y, z coordinates
    """
    x = r * np.cos(phi) * np.sin(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(theta)
    return np.array([x, y, z])


@njit
def cartesian_to_spherical(x, y, z):
    """Converts Cartesian coordinates to spherical

    Parameters
    ----------
    x : float64
        x coordinate
    y : float64
        y coordinate
    z : float64
        z coordinate

    Returns
    -------
    float64, float64, float64
        r, theta, phi coordinates
    """
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    return r, theta, phi


@njit
def doppler_gamma(direction_vector, position_vector, time):
    """Doppler shift for photons in 3D

    Parameters
    ----------
    direction_vector : numpy.ndarray
        Array of r, theta, phi vector (length 3)
    ejecta_velocity : float64
        Velocity of the ejecta

    Returns
    -------
    float64
        Doppler factor
    """
    velocity_vector = position_vector / time
    return 1 - (np.dot(direction_vector, velocity_vector) / C_CGS)


@njit
def angle_aberration_gamma(direction_vector, position_vector, time):
    """Angle aberration formula for photons in 3D

    Parameters
    ----------
    direction_vector : numpy.ndarray
        Array of r, theta, phi vector (length 3)
    ejecta_velocity : float64
        Velocity of the ejecta

    Returns
    -------
    numpy.ndarray
        New direction after aberration
    """
    velocity_vector = position_vector / time
    direction_dot_velocity = np.dot(direction_vector, velocity_vector)

    gamma = 1 / np.sqrt(1.0 - (np.dot(velocity_vector, velocity_vector) / (C_CGS ** 2)))

    factor_a = gamma * (1 - direction_dot_velocity / C_CGS)

    factor_b = (
        gamma - (gamma ** 2 * direction_dot_velocity / (gamma + 1) / C_CGS)
    ) / C_CGS

    output_vector = direction_vector - (velocity_vector * factor_b) / factor_a

    return np.array([1, output_vector[1], output_vector[2]])


@njit
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
    return energy / ELECTRON_MASS_ENERGY_KEV


@njit
def euler_rodrigues(theta, direction):
    """
    Calculates the Euler-Rodrigues rotation matrix

    Parameters
    ----------
    theta : float
        angle of rotation in radians
    direction : One dimensional Numpy array, dtype float
        x, y, z direction vector

    Returns
    -------
    rotation matrix : Two dimensional Numpy array, dtype float

    """
    a = np.cos(theta / 2)
    b = direction[0] * np.sin(theta / 2)
    c = direction[1] * np.sin(theta / 2)
    d = direction[2] * np.sin(theta / 2)

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


@njit
def solve_quadratic_equation(x, y, z, x_dir, y_dir, z_dir, radius):
    """
    Solves the quadratic equation for the distance to the shell boundary

    Parameters
    ----------
    x,y,z : float
    x_dir, y_dir, z_dir : float
    radius : float

    Returns
    -------
    solution_1 : float
    solution_2 : float

    """
    b = 2.0 * (x * x_dir + y * y_dir + z * z_dir)
    c = -(radius ** 2) + x ** 2 + y ** 2 + z ** 2
    root = b ** 2 - 4 * c
    solution_1 = -np.inf
    solution_2 = -np.inf
    if root > 0.0:
        solution_1 = 0.5 * (-b + np.sqrt(root))
        solution_2 = 0.5 * (-b - np.sqrt(root))
    elif root == 0:
        solution_1 = -0.5 * b

    return solution_1, solution_2


@njit
def klein_nishina(energy, theta_C):
    """
    Calculates the Klein-Nishina equation

    https://en.wikipedia.org/wiki/Klein%E2%80%93Nishina_formula

    .. math::
        \frac{r_e}{2} [1 + \kappa (1 - \cos\theta_C)]^{-2} \left( 1 + \cos^2\theta_C + \frac{\kappa^2 (1 - \cos\theta_C)^2}{1 + \kappa(1 - \cos\theta_C)}\right)

    where :math:`\kappa = E / (m_e c^2)`

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
        R_ELECTRON_SQUARED
        / 2
        * (1.0 + kappa * (1.0 - np.cos(theta_C))) ** -2.0
        * (
            1.0
            + np.cos(theta_C) ** 2.0
            + (kappa ** 2.0 * (1.0 - np.cos(theta_C)) ** 2.0)
            / (1.0 + kappa * (1.0 - np.cos(theta_C)))
        )
    )


@njit
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
    theta_angles = np.linspace(0, np.pi, sample_resolution)

    theta_distribution = np.cumsum(klein_nishina(energy, theta_angles))
    norm_theta_distribution = theta_distribution / np.max(theta_distribution)

    return theta_angles, norm_theta_distribution


@njit
def get_random_theta_photon():
    """Get a random theta direction between 0 and pi
    Returns
    -------
    float
        Random theta direction
    """
    return np.arccos(1.0 - 2.0 * np.random.random())


@njit
def get_random_phi_photon():
    """Get a random phi direction between 0 and 2 * pi

    Returns
    -------
    float
        Random phi direction
    """
    return 2.0 * np.pi * np.random.random()


@njit
def get_random_theta_photon_array(n):
    """Get a random theta direction between 0 and pi
    Returns
    -------
    float
        Random theta direction
    """
    return np.arccos(1.0 - 2.0 * np.random.random(n))


@njit
def get_random_phi_photon_array(n):
    """Get a random phi direction between 0 and 2 * pi

    Returns
    -------
    float
        Random phi direction
    """
    return 2.0 * np.pi * np.random.random(n)


def convert_half_life_to_astropy_units(half_life_string):
    """Converts input half-life to use astropy units

    Parameters
    ----------
    half_life_string : string
        Half-life as a string

    Returns
    -------
    astropy Quantity
        Half-life in seconds
    """
    value, unit = half_life_string.split(" ")
    try:
        nominal_value, magnitude = value.split("×")
        base, exponent = magnitude.split("+")
        nominal_value = float(nominal_value) * float(base) ** float(exponent)
    except ValueError:
        nominal_value = float(value)
    if unit == "y":
        unit = "yr"
    half_life_with_unit = nominal_value * u.Unit(unit)
    return half_life_with_unit.to(u.s)


@njit
def normalize_vector(vector):
    """
    Normalizes a vector in cartesian coordinates

    Parameters
    ----------
    vector : One-dimensional Numpy Array, dtype float
        Input vector

    Returns
    -------
    One-dimensional Numpy Array, dtype float
        Normalized vector
    """
    return vector / np.linalg.norm(vector)


@njit
def get_perpendicular_vector(original_direction):
    """
    Computes a vector which is perpendicular to the input vector

    Parameters
    ----------
    original_direction : SphericalVector object

    Returns
    -------
    numpy.ndarray
        Perpendicular vector to the input
    """
    # draw random angles
    theta = get_random_theta_photon()
    phi = get_random_phi_photon()
    # transform random angles to cartesian coordinates
    random_vector = spherical_to_cartesian(1, theta, phi)
    perpendicular_vector = np.cross(original_direction, random_vector)
    return normalize_vector(perpendicular_vector)


@njit
def get_index(value, array):
    """Get the index that places a value
    at array[i] < array <= vec[i+1]

    Parameters
    ----------
    value : float
        Value to locate
    array : array
        Array to search

    Returns
    -------
    int
        Index
    """    
    if value <= array[0]:
        return 0
    elif value > array[-1]:
        return len(array) - 1

    i = 0
    while value > array[i+1]:
        i += 1

    return i