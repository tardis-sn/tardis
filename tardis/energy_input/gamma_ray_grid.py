import numpy as np
from tardis.energy_input.util import solve_quadratic_equation
from astropy.coordinates import cartesian_to_spherical


def calculate_distance_radial(gxpacket, r_inner, r_outer):
    """
    Calculates 3D distance to shell from gamma ray position

    Parameters
    ----------
    gxpacket : GXPacket object
    r_inner : dtype float
    r_outer : dtype float

    Returns
    -------
    distance : dtype float

    """
    # determine cartesian location coordinates of gamma-ray object
    x, y, z = gxpacket.location.get_cartesian_coords
    # determine cartesian direction coordinates of gamma-ray object
    x_dir, y_dir, z_dir = gxpacket.direction.get_cartesian_coords
    # solve the quadratic distance equation for the inner and
    # outer shell boundaries
    inner_1, inner_2 = solve_quadratic_equation(
        x, y, z, x_dir, y_dir, z_dir, r_inner
    )
    outer_1, outer_2 = solve_quadratic_equation(
        x, y, z, x_dir, y_dir, z_dir, r_outer
    )
    distances = [inner_1, inner_2, outer_1, outer_2]
    # the correct distance is the shortest positive distance
    distance = min(i for i in distances if i > 0.0)

    return distance


def distance_trace(
    gxpacket,
    inner_radii,
    outer_radii,
    total_opacity,
    ejecta_epoch,
):
    """
    Traces distance traveled by gamma ray and finds distance to
    next interaction and boundary

    Parameters
    ----------
    gxpacket : GXPacket object
    inner_radii : One dimensional Numpy array, dtype float
    outer_radii : One dimensional Numpy array, dtype float
    total_opacity : dtype float
    ejecta_epoch : dtype float

    Returns
    -------
    distance_interaction : dtype float
    distance_boundary : dtype float

    """
    if gxpacket.shell < len(inner_radii):
        distance_boundary = calculate_distance_radial(
            gxpacket,
            inner_radii[gxpacket.shell],
            outer_radii[gxpacket.shell],
        )
    else:
        distance_boundary = 0.0

    distance_interaction = gxpacket.tau / total_opacity / ejecta_epoch

    return distance_interaction, distance_boundary


def move_gamma_ray(gxpacket, distance):
    """
    Moves gamma ray a distance along its direction vector

    Parameters
    ----------
    gxpacket : GXPacket object
    distance : dtype float

    Returns
    -------
    gxpacket : GXPacket object

    """
    x_old, y_old, z_old = gxpacket.location.get_cartesian_coords
    x_dir, y_dir, z_dir = gxpacket.direction.get_cartesian_coords
    # overshoot by 1e-7 * distance to shell boundary
    # so that the gamma-ray is comfortably in the next shell
    x_new = x_old + distance * (1 + 1e-7) * x_dir
    y_new = y_old + distance * (1 + 1e-7) * y_dir
    z_new = z_old + distance * (1 + 1e-7) * z_dir

    r, theta, phi = cartesian_to_spherical(x_new, y_new, z_new)
    gxpacket.location.r = r.value
    gxpacket.location.theta = theta.value + 0.5 * np.pi
    gxpacket.location.phi = phi.value
    return gxpacket


def density_sampler(radii, mass_ratio):
    """
    Randomly samples the

    Parameters
    ----------
    radii : GammaRay object
    mass_ratio : dtype float

    Returns
    -------
    radius : dtype float
    index : dtype int

    """
    z = np.random.random()

    mass_ratio_sorted_indices = np.argsort(mass_ratio)
    index = mass_ratio_sorted_indices[
        np.searchsorted(mass_ratio, z, sorter=mass_ratio_sorted_indices)
    ]

    return radii[index], index


def mass_distribution(
    radial_grid_size, inner_radii, outer_radii, density_profile
):
    """Calculates the distribution of mass in the shells
    based on a density profile

    Parameters
    ----------
    radial_grid_size : int
        Number of radial grid cells
    inner_radii : One-dimensional Numpy Array, dtype float
        Inner radii of shells
    outer_radii : One-dimensional Numpy Array, dtype float
        Outer radii of shells
    density_profile : One-dimensional Numpy Array, dtype float
        Density of shells

    Returns
    -------
    One-dimensional Numpy Array, dtype float
        Normalized array of mass in each shell
    """
    mass = np.zeros(radial_grid_size)

    i = 0
    while i < radial_grid_size:
        if i == 0:
            mass[i] = (
                4.0
                / 3.0
                * np.pi
                * density_profile[i]
                * (outer_radii[i] - inner_radii[i]) ** 3.0
            )
        else:
            mass[i] = (
                4.0
                / 3.0
                * np.pi
                * density_profile[i]
                * (outer_radii[i] ** 3.0 - outer_radii[i - 1] ** 3.0)
            )

            mass[i] += mass[i - 1]

        i += 1

    return mass / np.max(mass)


def get_shell(radius, outer_radii):
    """Returns the shell index at a given radius

    Parameters
    ----------
    radius : float
        Radius of interest
    outer_radii : One-dimensional Numpy Array, dtype float
        Outer radii of shells

    Returns
    -------
    int
        Shell index corresponding to radius
    """
    shell_inner = np.searchsorted(outer_radii, radius, side="left")

    return shell_inner
