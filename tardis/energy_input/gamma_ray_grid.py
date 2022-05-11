import numpy as np
import re
from nuclear.io.nndc import get_decay_radiation_database, store_decay_radiation
import pandas as pd
import astropy.units as u
from numba import njit

from tardis.energy_input.util import (
    doppler_gamma,
    solve_quadratic_equation,
    convert_half_life_to_astropy_units,
    cartesian_to_spherical,
    C_CGS,
)
from tardis.util.base import (
    atomic_number2element_symbol,
)
import tardis.constants as const


@njit
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

    # solve the quadratic distance equation for the inner and
    # outer shell boundaries
    inner_1, inner_2 = solve_quadratic_equation(
        photon.location, photon.direction, r_inner
    )
    outer_1, outer_2 = solve_quadratic_equation(
        photon.location, photon.direction, r_outer
    )
   
    final_position_inner_1 = photon.location + photon.direction*inner_1
    final_position_inner_2 = photon.location + photon.direction*inner_2
    final_position_outer_1 = photon.location + photon.direction*outer_1
    final_position_outer_2 = photon.location + photon.direction*outer_2

    if np.dot(final_position_inner_1, photon.direction) > 0:
        inner_1 = -1
    if np.dot(final_position_inner_2, photon.direction) > 0:
        inner_2 = -1
    if np.dot(final_position_outer_1, photon.direction) < 0:
        outer_1 = -1
    if np.dot(final_position_outer_2, photon.direction) < 0:
        outer_2 = -1

    distances = np.array([inner_1, inner_2, outer_1, outer_2])

    # the correct distance is the shortest positive distance
    distance_list = [i for i in distances if i > 0]

    if not distance_list:
        print(photon.get_location_r() - r_inner)
        print(photon.get_location_r() - r_outer)
        print(photon.get_location_r())
        print(photon.location, photon.direction, r_inner, r_outer)
        print(distances)
        print(photon.shell)
        raise ValueError("No root found for distance calculation!")

    shortest = min(distance_list)
    shell_change = 1

    if shortest == (inner_1 or inner_2):
        shell_change = -1

    return shortest, shell_change


@njit
def distance_trace(
    photon,
    inner_velocity,
    outer_velocity,
    total_opacity,
    current_time,
    next_time
):
    """
    Traces distance traveled by gamma ray and finds distance to
    next interaction and boundary

    Parameters
    ----------
    photon : GXPhoton object
    inner_velocity : One dimensional Numpy array, dtype float
    outer_velocity : One dimensional Numpy array, dtype float
    total_opacity : float
    current_time : float
    next_time : float

    Returns
    -------
    distance_interaction : float
    distance_boundary : float
    distance_time : float
    shell_change : int
    """
    distance_boundary, shell_change = calculate_distance_radial(
        photon,
        inner_velocity[photon.shell] * current_time,
        outer_velocity[photon.shell] * current_time,
    )

    distance_interaction = photon.tau / total_opacity
    distance_time = (next_time - photon.time_current) * C_CGS
    return distance_interaction, distance_boundary, distance_time, shell_change


@njit
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
    x_old, y_old, z_old = photon.get_location_cartesian_coords()
    x_dir, y_dir, z_dir = photon.get_direction_cartesian_coords()

    y_new = y_old + distance * y_dir
    z_new = z_old + distance * z_dir
    x_new = x_old + distance * x_dir

    r, theta, phi = cartesian_to_spherical(x_new, y_new, z_new)
    photon.location_r = r
    photon.location_theta = theta
    photon.location_phi = phi

    return photon


@njit
def move_packet(packet, distance):
    """
    Moves packet a distance along its direction vector

    Parameters
    ----------
    packet : GXPacket object
    distance : float

    Returns
    -------
    packet : GXPacket object

    """
    location_old = packet.location
    direction = packet.direction

    location_new = location_old + distance * direction

    packet.location = location_new

    doppler_factor = doppler_gamma(
        packet.direction, packet.location, packet.time_current
    )

    packet.nu_cmf = packet.nu_rf * doppler_factor
    packet.energy_cmf = packet.energy_rf * doppler_factor

    return packet


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
    number_of_decays : int64
        Total number of simulation decays

    Returns
    -------
    pandas.DataFrame
        Photons required per shell
    pandas.DataFrame
        Database of decay radiation
    pandas.DataFrame
        Activity per shell
    pandas.DataFrame
        Activity scaled by decays per shell
    """

    abundance_dict = {}
    nuclide_mass_dict = {}
    for isotope_string, row in isotope_abundance.iterrows():
        if isotope_string == "Fe56":
            continue
        store_decay_radiation(isotope_string, force_update=False)
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
    scaled_activity_df = activity_df.copy()

    for column in decays_per_shell_df:
        scaled_decays_per_shell = (
            decays_per_shell_df[column]
            / abundance_norm_total_activity
            * number_of_decays
        )
        decays_per_shell_df[column] = round(scaled_decays_per_shell).astype(int)
        scaled_activity_df[column] /= scaled_decays_per_shell

    return (
        decays_per_shell_df,
        decay_rad_db,
        activity_per_shell,
        scaled_activity_df,
        activity_df,
    )


def compute_required_packets_per_shell(
    shell_masses,
    isotope_abundance,
    number_of_packets,
):
    """Computes the number of packets required per shell
    that sum to the total number of requested packets.
    Also stores/updates decay radiation in an HDF file.

    Parameters
    ----------
    shell_masses : ndarray
        Array of shell masses
    isotope_abundance : pandas DataFrame
        Abundances of isotopes
    number_of_decays : int64
        Total number of requested packets

    Returns
    -------
    pandas.DataFrame
        Packets required per shell
    pandas.DataFrame
        Database of decay radiation
    pandas.DataFrame
        Activity per shell
    pandas.DataFrame
        Activity scaled by decays per shell
    """

    for isotope_string, row in isotope_abundance.iterrows():
        if isotope_string == "Fe56":
            continue
        store_decay_radiation(isotope_string, force_update=False)

    unstable_isotope_abundance = isotope_abundance.drop("Fe56", inplace=False)

    mass_fraction_df = unstable_isotope_abundance.copy()
    nuclide_mass_df = unstable_isotope_abundance.copy() * shell_masses

    decay_rad_db, meta = get_decay_radiation_database()

    activity_df = mass_fraction_df.T.copy()
    for column in activity_df:
        isotope_meta = meta.loc[column]
        half_life = isotope_meta.loc[
            isotope_meta["key"] == "Parent T1/2 value"
        ]["value"].values[0]
        half_life = convert_half_life_to_astropy_units(half_life)
        decay_constant = np.log(2) / half_life
        atomic_mass = float(re.findall("\d+", column)[0]) * u.u.to(
            u.g / u.mol, equivalencies=u.molar_mass_amu()
        )
        number_of_nuclides = (
            nuclide_mass_df.T[column] / atomic_mass
        ) * const.N_A

        activity_df[column] = decay_constant * number_of_nuclides

    total_activity = activity_df.to_numpy().sum()
    decays_per_shell_df = activity_df.copy()
    scaled_activity_df = activity_df.copy()

    for column in decays_per_shell_df:
        scaled_decays_per_shell = (
            activity_df[column] / total_activity * number_of_packets
        )
        decays_per_shell_df[column] = round(scaled_decays_per_shell).astype(int)
        scaled_activity_df[column] /= decays_per_shell_df[
            column
        ]  # scaled_decays_per_shell

    return (
        decays_per_shell_df,
        scaled_activity_df,
        activity_df,
    )

def mass_fraction_packets_per_shell(
    isotope_abundance,
    number_of_packets,
):
    """Calculates packets per shell by mass fraction

    Parameters
    ----------
    isotope_abundance : DataFrame
        Isotope abundance dataframe
    number_of_packets : int
        Number of packets

    Returns
    -------
    DataFrame
        Packets per shell by mass fraction
    """    
    if "Fe56" in isotope_abundance.columns:
        mass_fraction_df = isotope_abundance.drop(columns="Fe56", inplace=False)
    else:
        mass_fraction_df = isotope_abundance

    total_mass_fraction_of_nuclides = mass_fraction_df.sum().sum()

    for column in mass_fraction_df:
        mass_fraction_df[column] = round((mass_fraction_df[column] / total_mass_fraction_of_nuclides) * number_of_packets).astype(int)

    return mass_fraction_df

def activity_per_shell(
    isotope_masses
):  
    """Calculate isotope activity per shell based on mass

    Parameters
    ----------
    isotope_masses : DataFrame
        DataFrame of simulation isotope masses per shell

    Returns
    -------
    DataFrame
        Activity per shell
    DataFrame
        Decay radiation database
    DataFrame
        Metadata for the decay radiation database
    """    
    for isotope_string in isotope_masses:
        if isotope_string == "Fe56":
            isotope_masses.drop(columns="Fe56", inplace=True)
            continue
        #store_decay_radiation(isotope_string, force_update=False)

    nuclide_mass_df = isotope_masses.copy()

    decay_rad_db, meta = get_decay_radiation_database()

    activity_df = nuclide_mass_df.copy()
    for column in activity_df:
        isotope_meta = meta.loc[column]
        half_life = isotope_meta.loc[
            isotope_meta["key"] == "Parent T1/2 value"
        ]["value"].values[0]
        half_life = convert_half_life_to_astropy_units(half_life)
        decay_constant = np.log(2) / half_life
        atomic_mass = float(re.findall("\d+", column)[0]) * u.u.to(
            u.g / u.mol, equivalencies=u.molar_mass_amu()
        )
        number_of_nuclides = (
            nuclide_mass_df[column] / atomic_mass
        ) * const.N_A

        activity_df[column] = decay_constant * number_of_nuclides

    return activity_df, decay_rad_db, meta



def get_decay_database(
    isotope_abundance,
):
    """Gets the decay radiation database for a set
    of isotopes

    Parameters
    ----------
    isotope_abundance : DataFrame
        DataFrame of simulation isotope masses per shell

    Returns
    -------
    DataFrame
        Decay radiation database
    DataFrame
        Metadata for the decay radiation database
    """    
    for column in isotope_abundance:
        if column == "Fe56":
            continue
        store_decay_radiation(column, force_update=False)

    decay_rad_db, meta = get_decay_radiation_database()

    return decay_rad_db, meta


def get_tau(meta, isotope_string):
    """Calculate the mean lifetime of an isotope

    Parameters
    ----------
    meta : DataFrame
        Isotope metadata
    isotope_string : str
        Isotope of interest

    Returns
    -------
    float
        Mean lifetime of isotope
    """    
    isotope_meta = meta.loc[isotope_string]
    half_life = isotope_meta.loc[
        isotope_meta["key"] == "Parent T1/2 value"
    ]["value"].values[0]
    half_life = convert_half_life_to_astropy_units(half_life)
    return half_life / np.log(2)

def get_isotope_string(atom_number, atom_mass):
    """Get the isotope string in the format e.g. Ni56

    Parameters
    ----------
    atom_number : int
        Atomic number
    atom_mass : int
        Atomic mass

    Returns
    -------
    str
        Isotope string in the format e.g. Ni56
    """    
    return atomic_number2element_symbol(atom_number) + str(atom_mass)

def read_artis_lines(isotope):
    """Reads lines of ARTIS format

    Parameters
    ----------
    isotope : string
        Isotope to read e.g. Ni56

    Returns
    -------
    pd.DataFrame
        Energies and intensities of the isotope lines
    """    
    return pd.read_csv(
        "~/Downloads/tardisnuclear/" + isotope + ".txt",
        names=["energy", "intensity"],
        sep="  ",
        index_col=False,
    )
