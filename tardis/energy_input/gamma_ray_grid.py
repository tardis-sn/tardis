import numpy as np
from tardis.energy_input.util import solve_quadratic_equation
from astropy.coordinates import cartesian_to_spherical


def calculate_distance_radial(gamma_ray, r_inner, r_outer):
    """
    Calculates 3D distance to shell from gamma ray position

    Parameters
    ----------
    gamma_ray : GammaRay object
    r_inner : dtype float
    r_outer : dtype float

    Returns
    -------
    distance : dtype float

    """
    # determine cartesian location coordinates of gamma-ray object
    x, y, z = gamma_ray.location.get_cartesian_coords
    # determine cartesian direction coordinates of gamma-ray object
    x_dir, y_dir, z_dir = gamma_ray.direction.get_cartesian_coords
    # solve the quadratic distance equation for the inner and
    # outer shell boundaries
    inner_1, inner_2 = solve_quadratic_equation(
        x, y, z, x_dir, y_dir, z_dir, r_inner
    )
    outer_1, outer_2 = solve_quadratic_equation(
        x, y, z, x_dir, y_dir, z_dir, r_outer
    )
    distances = [inner_1, inner_2, outer_1, outer_2]
    # print(distances)
    # the correct distance is the shortest positive distance
    distance = min(i for i in distances if i > 0.0)
    # print("selected distance ", distance)

    return distance


def distance_trace(
    gamma_ray,
    inner_radii,
    outer_radii,
    total_opacity,
    distance_moved,
    ejecta_epoch,
):
    """
    Traces distance traveled by gamma ray and finds distance to
    next interaction and boundary

    Parameters
    ----------
    gamma_ray : GammaRay object
    radii : One dimensional Numpy array, dtype float
    total_opacity : dtype float
    distance_moved : dtype float

    Returns
    -------
    distance_interaction : dtype float
    distance_boundary : dtype float
     : dtype bool

    """
    if gamma_ray.shell < len(inner_radii):
        distance_boundary = calculate_distance_radial(
            gamma_ray,
            inner_radii[gamma_ray.shell],
            outer_radii[gamma_ray.shell],
        )
    else:
        distance_boundary = 0.0
    # print("current tau: ", gamma_ray.tau)
    distance_interaction = gamma_ray.tau / total_opacity / ejecta_epoch
    # print("distance interaction: ", distance_interaction)
    # print("distance boundary: ", distance_boundary)
    if distance_interaction < distance_boundary:
        gamma_ray.tau -= total_opacity * distance_interaction * ejecta_epoch
        # print("remaining tau: ", gamma_ray.tau)
        return distance_interaction, distance_boundary, True
    else:
        gamma_ray.tau -= total_opacity * distance_boundary * ejecta_epoch
        # print("remaining tau: ", gamma_ray.tau)
        return distance_interaction, distance_boundary, False


def move_gamma_ray(gamma_ray, distance):
    """
    Moves gamma ray a distance along its direction vector

    Parameters
    ----------
    gamma_ray : GammaRay object
    distance : dtype float

    Returns
    -------
    gamma_ray : GammaRay object

    """
    x_old, y_old, z_old = gamma_ray.location.get_cartesian_coords
    x_dir, y_dir, z_dir = gamma_ray.direction.get_cartesian_coords
    # print("direction vector: ", x_dir, y_dir, z_dir)
    # overshoot by 1e-7 * distance to shell boundary
    # so that the gamma-ray is comfortably in the next shell
    x_new = x_old + distance * (1 + 1e-7) * x_dir
    y_new = y_old + distance * (1 + 1e-7) * y_dir
    z_new = z_old + distance * (1 + 1e-7) * z_dir
    # print("selected distance: ", distance)
    # print("velocity coordinate change: ")
    # print("x:\t", distance * (1 + 1e-7) * x_dir)
    # print("y:\t", distance * (1 + 1e-7) * y_dir)
    # print("z:\t", distance * (1 + 1e-7) * z_dir)
    # print(
    #     "moved from (",
    #     round(x_old, 5),
    #     round(y_old, 5),
    #     round(z_old, 5),
    #     ")\n\t to (",
    #     round(x_new, 5),
    #     round(y_new, 5),
    #     round(z_new, 5),
    #     ")",
    # )
    r, theta, phi = cartesian_to_spherical(x_new, y_new, z_new)
    gamma_ray.location.r = r.value
    gamma_ray.location.theta = theta.value + 0.5 * np.pi
    gamma_ray.location.phi = phi.value
    return gamma_ray


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
