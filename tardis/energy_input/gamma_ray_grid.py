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
import tardis.constants as const


def calculate_distance_radial(photon, r_inner, r_outer):
    """
    Calculates 3D distance to shell from gamma ray position

    Parameters
    ----------
    photon : GXPhoton object
    r_inner : float
    r_outer : float

    Returns
    -------
    distance : float

    """
    # TODO: Maybe only calculate distances that are strictly needed instead of all four by default?
    # determine cartesian location coordinates of gamma-ray object
    x, y, z = photon.location.cartesian_coords
    # determine cartesian direction coordinates of gamma-ray object
    x_dir, y_dir, z_dir = photon.direction.cartesian_coords
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
    photon,
    inner_radii,
    outer_radii,
    total_opacity,
    time_explosion,
):
    """
    Traces distance traveled by gamma ray and finds distance to
    next interaction and boundary

    Parameters
    ----------
    photon : GXPhoton object
    inner_radii : One dimensional Numpy array, dtype float
    outer_radii : One dimensional Numpy array, dtype float
    total_opacity : float
    time_explosion : float

    Returns
    -------
    distance_interaction : float
    distance_boundary : float

    """
    distance_boundary = calculate_distance_radial(
        photon,
        inner_radii[photon.shell],
        outer_radii[photon.shell],
    )

    distance_interaction = photon.tau / total_opacity / time_explosion
    return distance_interaction, distance_boundary


def move_photon(photon, distance):
    """
    Moves gamma ray a distance along its direction vector

    Parameters
    ----------
    photon : GXPhoton object
    distance : float

    Returns
    -------
    photon : GXPhoton object

    """
    x_old, y_old, z_old = photon.location.cartesian_coords
    x_dir, y_dir, z_dir = photon.direction.cartesian_coords

    y_new = y_old + distance * y_dir
    z_new = z_old + distance * z_dir
    x_new = x_old + distance * x_dir

    r, theta, phi = cartesian_to_spherical(x_new, y_new, z_new)
    photon.location.r = r.value
    # Plus 0.5 * pi to correct for astropy rotation frame
    photon.location.theta = theta.value + 0.5 * np.pi
    photon.location.phi = phi.value
    return photon


def compute_required_photons_per_shell(
    shell_masses,
    raw_isotope_abundance,
    number_of_photons,
):
    """Computes the number of photons required per shell
    that sum to the total number of requested photons.
    Also stores/updates decay radiation in an HDF file.

    Parameters
    ----------
    shell_masses : ndarray
        Array of shell masses
    raw_isotope_abundance : pandas DataFrame
        Abundances of isotopes
    number_of_photons : int64
        Total number of simulation photons

    Returns
    -------
    pandas DataFrame
        Photons required per shell
    pandas DataFrame
        Database of decay radiation
    """

    norm_shell_masses = shell_masses / np.sum(shell_masses)
    abundance_dict = {}
    nuclide_mass_dict = {}
    for (atom_number, atom_mass), row in raw_isotope_abundance.iterrows():
        isotope_string = atomic_number2element_symbol(atom_number) + str(
            atom_mass
        )
        store_decay_radiation(isotope_string, force_update=False)
        abundance_dict[isotope_string] = row * norm_shell_masses
        nuclide_mass_dict[isotope_string] = row * shell_masses

    abundance_df = pd.DataFrame.from_dict(abundance_dict)
    nuclide_mass_df = pd.DataFrame.from_dict(nuclide_mass_dict)

    decay_rad_db, meta = get_decay_radiation_database()

    activity_df = abundance_df.copy()
    decay_rate_per_shell_df = nuclide_mass_df.copy()
    for column in activity_df:
        isotope_meta = meta.loc[column]
        half_life = isotope_meta.loc[
            isotope_meta["key"] == "Parent T1/2 value"
        ]["value"].values[0]
        half_life = convert_half_life_to_astropy_units(half_life)
        atomic_mass = float(re.findall("\d+", column)[0])
        activity_factor = np.log(2) / atomic_mass / half_life
        activity_df[column] = activity_df[column] * activity_factor

        decay_rate_per_shell_df[column] = (
            (np.log(2) / half_life)
            * const.N_A
            * (nuclide_mass_df[column] / (atomic_mass ** 2.0))
        )

    total_activity = activity_df.to_numpy().sum()
    decay_rate_per_shell = decay_rate_per_shell_df.to_numpy().sum(axis=1)
    photon_per_shell_df = activity_df.copy()

    for column in photon_per_shell_df:
        photon_per_shell_df[column] = round(
            photon_per_shell_df[column] * number_of_photons / total_activity
        )
        photon_per_shell_df[column] = photon_per_shell_df[column].astype(int)

    return photon_per_shell_df, decay_rad_db, decay_rate_per_shell
