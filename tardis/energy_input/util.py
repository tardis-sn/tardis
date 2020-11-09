import astropy.units as u
import numpy as np

R_ELECTRON = 2.8179403227e-15

def kappa_calculation(energy):
    """
    Calculates kappa for various other calculations
    i.e. energy normalized to 511 KeV

    Parameters
    ----------
    energy : dtype float

    Returns
    -------
    kappa : dtype float

    """
    k = energy / (511.0e3 * u.eV)
    return k.value

def euler_rodrigues(theta, direction):
    """
    Calculates the Euler-Rodrigues rotation matrix

    Parameters
    ----------
    theta : dtype float
    direction : SphericalVector object

    Returns
    -------
    rotation matrix : Two dimensional Numpy array, dtype float

    """
    a = np.cos(theta / 2)
    b = direction.x * np.sin(theta / 2)
    c = direction.y * np.sin(theta / 2)
    d = direction.z * np.sin(theta / 2)
    
    er11 = a ** 2. + b ** 2. - c ** 2. - d ** 2.
    er12 = 2. * (b * c - a * d)
    er13 = 2. * (b * d + a * d)
    
    er21 = 2. * (b * c + a * d)
    er22 = a ** 2. + c ** 2. - b ** 2. - d ** 2.
    er23 = 2. * (c * d - a * b)
    
    er31 = 2. * (b * d - a * c)
    er32 = 2. * (c * d + a * b)
    er33 = a ** 2. + d ** 2. - b ** 2. - c ** 2.
    
    return np.array([[er11, er12, er13],
                     [er21, er22, er23],
                     [er31, er32, er33]])

def quadratic(b, c):
    """
    Solves the reduced quadratic equation

    Parameters
    ----------
    b : dtype float
    c : dtype float

    Returns
    -------
    x1, x2 : dtype float

    """
    delta = b ** 2 - 4.0 * c
    if delta > 0:
        delta = np.sqrt(delta)
        delta = delta * np.sign(b)
        q = -0.5 * (b + delta)
        x1 = q
        x2 = c / q
    elif delta < 0:
        x1 = -np.inf
        x2 = -np.inf
    else:
        x1 = -2.0 * c / b
        x2 = -np.inf
    return x1, x2

def klein_nishina(energy, theta_C):
    """
    Calculates the Klein-Nishina equation

    Parameters
    ----------
    energy : dtype float
    theta_C : dtype float

    Returns
    -------
     : dtype float

    """
    kappa = kappa_calculation(energy) 
    return R_ELECTRON / 2 * \
            (1. + kappa * (1. - np.cos(theta_C))) ** -2. * \
            (1. + np.cos(theta_C) ** 2. + \
            (kappa ** 2. * (1. - np.cos(theta_C)) ** 2.) / (1. + kappa * (1. - np.cos(theta_C))))

def compton_theta_distribution(energy, sample_resolution = 100):
    """
    Calculates the distribution of theta angles 
    for Compton Scattering

    Parameters
    ----------
    energy : dtype float
    sample_resolution : dtype int

    Returns
    -------
    theta_angles : One dimensional Numpy array, dtype float
    norm_theta_distribution : One dimensional Numpy array, dtype float

    """
    dtheta = np.pi / sample_resolution
    
    theta_distribution = np.zeros(sample_resolution)
    theta_angles = np.ones(sample_resolution) * np.pi
    
    theta = 0.
    for i in range(sample_resolution - 1):
        theta_distribution[i + 1] = \
                theta_distribution[i] + klein_nishina(energy, theta)
        theta_angles[i] = theta
        theta += dtheta

    norm_theta_distribution = theta_distribution / np.max(theta_distribution)
    
    return theta_angles, norm_theta_distribution
