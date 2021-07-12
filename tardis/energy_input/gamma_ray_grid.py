import numpy as np
from astropy.coordinates import cartesian_to_spherical
import re
from nuclear.io.nndc import get_decay_radiation_database, store_decay_radiation
import pandas as pd

from tardis.energy_input.util import (
    solve_quadratic_equation,
    convert_half_life_to_astropy_units,
)
from tardis.util.base import atomic_number2element_symbol


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
    # TODO: Maybe only calculate distances that are strictly needed instead of all four by default?
    # determine cartesian location coordinates of gamma-ray object
    x, y, z = gxpacket.location.cartesian_coords
    # determine cartesian direction coordinates of gamma-ray object
    x_dir, y_dir, z_dir = gxpacket.direction.cartesian_coords
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
    distance_boundary = calculate_distance_radial(
        gxpacket,
        inner_radii[gxpacket.shell],
        outer_radii[gxpacket.shell],
    )

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
    x_old, y_old, z_old = gxpacket.location.cartesian_coords
    x_dir, y_dir, z_dir = gxpacket.direction.cartesian_coords

    y_new = y_old + distance * y_dir
    z_new = z_old + distance * z_dir
    x_new = x_old + distance * x_dir

    r, theta, phi = cartesian_to_spherical(x_new, y_new, z_new)
    gxpacket.location.r = r.value
    gxpacket.location.theta = theta.value + 0.5 * np.pi
    gxpacket.location.phi = phi.value
    return gxpacket


def density_sampler(radii, mass_ratio):
    """
    Randomly samples the density of the model shells

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


def compute_required_packets_per_shell(
    shell_masses,
    raw_isotope_abundance,
    number_of_packets,
):
    """Computes the number of packets required per shell
    that sum to the total number of requested packets.
    Also stores/updates decay radiation in an HDF file.

    Parameters
    ----------
    shell_masses : ndarray
        Array of shell masses
    raw_isotope_abundance : pandas DataFrame
        Abundances of isotopes
    number_of_packets : int64
        Total number of simulation packets

    Returns
    -------
    pandas DataFrame
        Packets required per shell
    pandas DataFrame
        Database of decay radiation
    """

    shell_masses = shell_masses / np.sum(shell_masses)
    abundance_dict = {}
    for index, row in raw_isotope_abundance.iterrows():
        isotope_string = atomic_number2element_symbol(index[0]) + str(index[1])
        store_decay_radiation(isotope_string, force_update=False)
        abundance_dict[isotope_string] = row * shell_masses

    abundance_df = pd.DataFrame.from_dict(abundance_dict)

    decay_rad_db, meta = get_decay_radiation_database()

    activity_df = abundance_df.copy()
    for column in activity_df:
        isotope_meta = meta.loc[column]
        half_life = isotope_meta.loc[
            isotope_meta["key"] == "Parent T1/2 value"
        ]["value"].values[0]
        half_life = convert_half_life_to_astropy_units(half_life)
        atomic_mass = float(re.findall("\d+", column)[0])
        activity_factor = np.log(2) / atomic_mass / half_life
        activity_df[column] = activity_df[column] * activity_factor

    total_activity = activity_df.to_numpy().sum()
    packet_per_shell_df = activity_df.copy()

    for column in packet_per_shell_df:
        packet_per_shell_df[column] = round(
            packet_per_shell_df[column] * number_of_packets / total_activity
        )
        packet_per_shell_df[column] = packet_per_shell_df[column].astype(int)

    return packet_per_shell_df, decay_rad_db
