import numpy as np
from astropy.coordinates import cartesian_to_spherical
import re
from nuclear.io.nndc import get_decay_radiation_database, store_decay_radiation
import pandas as pd
import astropy.units as u

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
    photon.location.r = r
    # Plus 0.5 * pi to correct for astropy rotation frame
    photon.location.theta = theta.value + 0.5 * np.pi
    photon.location.phi = phi.value
    return photon


def compute_required_photons_per_shell(
    shell_masses,
    isotope_abundance,
    number_of_decays,
):
    """Computes the number of photons required per shell
    that sum to the total number of requested photons.
    Also stores/updates decay radiation in an HDF file.

    Parameters
    ----------
    shell_masses : ndarray
        Array of shell masses
    isotope_abundance : pandas DataFrame
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

    abundance_dict = {}
    nuclide_mass_dict = {}
    for isotope_string, row in isotope_abundance.iterrows():
        if isotope_string == "Fe56":
            continue
        store_decay_radiation(isotope_string, force_update=True)
        abundance_dict[isotope_string] = row
        nuclide_mass_dict[isotope_string] = row * shell_masses

    abundance_df = pd.DataFrame.from_dict(abundance_dict)
    nuclide_mass_df = pd.DataFrame.from_dict(nuclide_mass_dict)

    decay_rad_db, meta = get_decay_radiation_database()

    abundance_norm_activity_df = abundance_df.copy()
    activity_df = nuclide_mass_df.copy()
    for column in abundance_norm_activity_df:
        isotope_meta = meta.loc[column]
        half_life = isotope_meta.loc[
            isotope_meta["key"] == "Parent T1/2 value"
        ]["value"].values[0]
        half_life = convert_half_life_to_astropy_units(half_life)
        decay_constant = np.log(2) / half_life
        atomic_mass = float(re.findall("\d+", column)[0]) * u.u.to(
            u.g / u.mol, equivalencies=u.molar_mass_amu()
        )
        number_of_nuclides = (nuclide_mass_df[column] / atomic_mass) * const.N_A

        abundance_norm_activity_df[column] *= decay_constant
        activity_df[column] = decay_constant * number_of_nuclides

    abundance_norm_total_activity = abundance_norm_activity_df.to_numpy().sum()
    activity_per_shell = activity_df.to_numpy().sum(axis=1)
    decays_per_shell_df = abundance_norm_activity_df.copy()

    for column in decays_per_shell_df:
        decays_per_shell_df[column] = round(
            decays_per_shell_df[column]
            / abundance_norm_total_activity
            * number_of_decays
        )
        decays_per_shell_df[column] = decays_per_shell_df[column].astype(int)

    return decays_per_shell_df, decay_rad_db, activity_per_shell
