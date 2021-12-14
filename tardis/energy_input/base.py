import numpy as np
import copy
from tqdm.auto import tqdm
import pandas as pd
import astropy.units as u

from tardis.energy_input.util import SphericalVector
from tardis.energy_input.GXPhoton import GXPhoton, GXPhotonStatus
from tardis.energy_input.gamma_ray_grid import (
    distance_trace,
    move_photon,
    compute_required_photons_per_shell,
)
from tardis.energy_input.energy_source import (
    setup_input_energy,
    sample_energy_distribution,
    intensity_ratio,
    decay_nuclides,
)
from tardis.energy_input.calculate_opacity import (
    compton_opacity_calculation,
    photoabsorption_opacity_calculation,
    pair_creation_opacity_calculation,
)
from tardis.energy_input.gamma_ray_interactions import (
    scatter_type,
    compton_scatter,
    pair_creation,
    get_compton_angle,
)
from tardis.energy_input.util import (
    get_random_theta_photon,
    get_random_phi_photon,
    cartesian_to_spherical,
    ELECTRON_MASS_ENERGY_KEV,
    BOUNDARY_THRESHOLD,
    KEV2ERG,
    C_CGS,
)
from tardis import constants as const
from tardis.util.base import (
    atomic_number2element_symbol,
)

# Energy: keV, exported as eV for SF solver
# distance: cm
# mass: g
# time: s


def initialize_photons(
    number_of_shells,
    decays_per_shell,
    ejecta_volume,
    shell_masses,
    inner_velocities,
    outer_velocities,
    decay_rad_db,
    scaled_activity_df,
):
    """Initializes photon properties
    and appends beta decay energy to output tables

    Parameters
    ----------
    number_of_shells : int
        Number of shells in model
    decays_per_shell : pandas.Dataframe
        Number of decays in a shell
    ejecta_volume : numpy.array
        Volume per shell
    shell_masses : numpy.array
        Mass per shell
    inner_velocities : numpy.array
        Shell inner velocities
    outer_velocities : numpy.array
        Shell outer velocities
    decay_rad_db : pandas.Dataframe
        Decay radiation database

    Returns
    -------
    list
        GXPhoton objects
    numpy array
        energy binned per shell
    list
        photon info
    """
    photons = []
    energy_df_rows = np.zeros(number_of_shells)
    energy_plot_df_rows = []

    cumulative_mass = np.cumsum(shell_masses)
    cumulative_mass /= np.max(cumulative_mass)

    cumulative_mass = np.insert(cumulative_mass, 0, 0)
    outer_velocities = np.insert(outer_velocities, 0, inner_velocities[0])

    scaled_decays_per_shell = decays_per_shell.copy()

    for column in scaled_decays_per_shell:
        subtable = decay_rad_db.loc[column]
        (
            gamma_ray_probability,
            positron_probability,
            scale_factor,
        ) = intensity_ratio(subtable, "'gamma_rays' or type=='x_rays'", "'e+'")
        energy_sorted, energy_cdf = setup_input_energy(
            subtable, "'gamma_rays' or type=='x_rays'"
        )
        positron_energy_sorted, positron_energy_cdf = setup_input_energy(
            subtable, "'e+'"
        )
        for shell in range(number_of_shells):
            scaled_decays_per_shell[column].iloc[shell] *= scale_factor
            requested_decays_per_shell = int(
                scaled_decays_per_shell[column].iloc[shell]
            )
            for _ in range(requested_decays_per_shell):
                # draw a random gamma-ray in shell
                primary_photon = GXPhoton(
                    location=0,
                    direction=0,
                    energy=1,
                    status=GXPhotonStatus.IN_PROCESS,
                    shell=0,
                    activity=scaled_activity_df[column].iloc[shell],
                )

                z = (
                    cumulative_mass[shell]
                    + (cumulative_mass[shell + 1] - cumulative_mass[shell])
                    * np.random.random()
                )

                initial_radius = np.interp(z, cumulative_mass, outer_velocities)

                location_theta = get_random_theta_photon()
                location_phi = get_random_phi_photon()
                primary_photon.location = SphericalVector(
                    initial_radius, location_theta, location_phi
                )
                direction_theta = get_random_theta_photon()
                direction_phi = get_random_phi_photon()
                primary_photon.direction = SphericalVector(
                    1.0, direction_theta, direction_phi
                )
                primary_photon.shell = shell
                if gamma_ray_probability < np.random.random():
                    # positron: sets gamma-ray energy to 511keV
                    primary_photon.energy = ELECTRON_MASS_ENERGY_KEV
                    photons.append(primary_photon)

                    # annihilation dumps energy into medium
                    energy_KeV = sample_energy_distribution(
                        positron_energy_sorted, positron_energy_cdf
                    )

                    # convert KeV to eV / cm^3
                    energy_df_rows[shell] += (
                        energy_KeV
                        * primary_photon.activity
                        * 1000
                        / ejecta_volume[shell]
                    )
                    energy_plot_df_rows.append(
                        [
                            -1,
                            energy_KeV,
                            initial_radius,
                            primary_photon.location.theta,
                            0.0,
                            -1,
                        ]
                    )

                    # annihilation produces second gamma-ray in opposite direction
                    (
                        primary_photon_x,
                        primary_photon_y,
                        primary_photon_z,
                    ) = primary_photon.direction.cartesian_coords
                    secondary_photon = copy.deepcopy(primary_photon)
                    secondary_photon_x = -primary_photon_x
                    secondary_photon_y = -primary_photon_y
                    secondary_photon_z = -primary_photon_z
                    (
                        secondary_photon_r,
                        secondary_photon_theta,
                        secondary_photon_phi,
                    ) = cartesian_to_spherical(
                        secondary_photon_x,
                        secondary_photon_y,
                        secondary_photon_z,
                    )
                    secondary_photon.direction.r = secondary_photon_r
                    secondary_photon.direction.theta = secondary_photon_theta
                    secondary_photon.direction.phi = secondary_photon_phi

                    photons.append(secondary_photon)
                else:
                    # Spawn a gamma ray emission with energy from gamma-ray list
                    primary_photon.energy = sample_energy_distribution(
                        energy_sorted, energy_cdf
                    )
                    photons.append(primary_photon)

    return photons, energy_df_rows, energy_plot_df_rows


def main_gamma_ray_loop(num_photons, model):
    """Main loop that determines the gamma ray propagation

    Parameters
    ----------
    num_photons : int
        Number of photons
    model : tardis.Radial1DModel
        The tardis model to calculate gamma ray propagation through

    Returns
    -------
    pandas.DataFrame
        Energy per mass per shell in units of eV/s/cm^-3
    pandas.DataFrame
        Columns:
        Photon index,
        Energy input per photon,
        radius of deposition,
        theta angle of deposition,
        time of deposition,
        type of deposition where:
            -1 = beta decay,
            0 = Compton scatter,
            1 = photoabsorption,
            2 = pair creation
    list
        Energy of escaping photons
    """
    escape_energy = []

    # Note the use of velocity as the radial coordinate
    # Enforce cgs
    outer_velocities = model.v_outer.to("cm/s").value
    inner_velocities = model.v_inner.to("cm/s").value
    ejecta_density = model.density.to("g/cm^3").value
    ejecta_volume = model.volume.to("cm^3").value
    time_explosion = model.time_explosion.to("s").value
    number_of_shells = model.no_of_shells
    raw_isotope_abundance = model.raw_isotope_abundance

    shell_masses = ejecta_volume * ejecta_density

    new_abundance_rows = []

    for index, mass in enumerate(shell_masses):
        isotope_abundance = raw_isotope_abundance[index]
        isotope_dict = {}
        for (
            atom_number,
            atom_mass,
        ), abundance in isotope_abundance.iteritems():
            isotope_string = atomic_number2element_symbol(atom_number) + str(
                atom_mass
            )
            isotope_dict[isotope_string] = abundance
        new_abundance_rows.append(
            decay_nuclides(
                mass * u.g.to("M_sun"),
                isotope_dict,
                np.array([time_explosion * u.s.to("day")]),
            )
        )

    new_abundances = pd.concat(new_abundance_rows).transpose()
    new_abundances.columns = raw_isotope_abundance.columns

    # scale abundances, needs to be done per-isotope in future
    new_abundances *= raw_isotope_abundance.values

    (
        decays_per_shell,
        decay_rad_db,
        decay_rate_per_shell,
        scaled_activity_df,
    ) = compute_required_photons_per_shell(
        shell_masses, new_abundances, num_photons
    )

    # Taking iron group to be elements 21-30
    # Used as part of the approximations for photoabsorption and pair creation
    # Dependent on atomic data
    iron_group_fraction_per_shell = model.abundance.loc[(21):(30)].sum(axis=0)

    (photons, energy_df_rows, energy_plot_df_rows,) = initialize_photons(
        number_of_shells,
        decays_per_shell,
        ejecta_volume,
        shell_masses,
        inner_velocities,
        outer_velocities,
        decay_rad_db,
        scaled_activity_df,
    )

    scaled_decay_rate_per_shell = (
        decay_rate_per_shell / decays_per_shell.to_numpy().sum(axis=1)
    )

    total_energy = np.zeros(number_of_shells)
    for p in photons:
        total_energy[p.shell] += p.energy * p.activity / ejecta_volume[p.shell]

    injected = (total_energy * KEV2ERG) + (energy_df_rows * KEV2ERG)

    print("=== Injected energy ===")
    print(injected)

    i = 0
    for photon in tqdm(photons):

        while photon.status == GXPhotonStatus.IN_PROCESS:

            compton_opacity = compton_opacity_calculation(
                photon.energy, ejecta_density[photon.shell]
            )
            photoabsorption_opacity = photoabsorption_opacity_calculation(
                photon.energy,
                ejecta_density[photon.shell],
                iron_group_fraction_per_shell[photon.shell],
            )
            pair_creation_opacity = pair_creation_opacity_calculation(
                photon.energy,
                ejecta_density[photon.shell],
                iron_group_fraction_per_shell[photon.shell],
            )
            total_opacity = (
                compton_opacity
                + photoabsorption_opacity
                + pair_creation_opacity
            )

            (distance_interaction, distance_boundary,) = distance_trace(
                photon,
                inner_velocities,
                outer_velocities,
                total_opacity,
                time_explosion,
            )
            if distance_interaction < distance_boundary:

                photon.tau = -np.log(np.random.random())

                photon.status = scatter_type(
                    compton_opacity,
                    photoabsorption_opacity,
                    total_opacity,
                )

                photon = move_photon(photon, distance_interaction)
                photon.shell = np.searchsorted(
                    outer_velocities, photon.location.r, side="left"
                )
                photon.time_current += (
                    distance_interaction / C_CGS * time_explosion
                )

                if photon.status == GXPhotonStatus.COMPTON_SCATTER:
                    (
                        compton_angle,
                        ejecta_energy_gained,
                        photon.energy,
                    ) = get_compton_angle(photon.energy)
                    (
                        photon.direction.theta,
                        photon.direction.phi,
                    ) = compton_scatter(photon, compton_angle)

                if photon.status == GXPhotonStatus.PAIR_CREATION:
                    ejecta_energy_gained = photon.energy - (
                        2.0 * ELECTRON_MASS_ENERGY_KEV
                    )
                    photon, backward_photon = pair_creation(photon)

                    # Add antiparallel photon on pair creation at end of list
                    photons.append(backward_photon)

                if photon.status == GXPhotonStatus.PHOTOABSORPTION:
                    ejecta_energy_gained = photon.energy

                # Save photons to dataframe rows
                # convert KeV to eV / s / cm^3
                energy_df_rows[photon.shell] += (
                    photon.activity
                    * ejecta_energy_gained
                    * 1000
                    / ejecta_volume[photon.shell]
                )
                energy_plot_df_rows.append(
                    [
                        i,
                        ejecta_energy_gained,
                        photon.location.r,
                        photon.location.theta,
                        photon.time_current,
                        int(photon.status),
                    ]
                )

                if photon.status == GXPhotonStatus.PHOTOABSORPTION:
                    # Photon destroyed, go to the next photon
                    break
                else:
                    photon.status = GXPhotonStatus.IN_PROCESS

            else:
                photon.tau -= total_opacity * distance_boundary * time_explosion
                # overshoot so that the gamma-ray is comfortably in the next shell
                photon = move_photon(
                    photon, distance_boundary * (1 + BOUNDARY_THRESHOLD)
                )
                photon.time_current += (
                    distance_boundary / C_CGS * time_explosion
                )
                photon.shell = np.searchsorted(
                    outer_velocities, photon.location.r, side="left"
                )

            if photon.shell > len(ejecta_density) - 1:
                escape_energy.append(photon.energy)
                photon.status = GXPhotonStatus.END
            elif photon.shell < 0:
                photon.energy = 0.0
                photon.status = GXPhotonStatus.END

        i += 1

    # DataFrame of energy information
    energy_plot_df = pd.DataFrame(
        data=energy_plot_df_rows,
        columns=[
            "photon_index",
            "energy_input",
            "energy_input_r",
            "energy_input_theta",
            "energy_input_time",
            "energy_input_type",
        ],
    )

    # Energy is eV/s/cm^-3
    energy_df = pd.DataFrame(data=energy_df_rows, columns=["energy"])

    return (
        energy_df,
        energy_plot_df,
        escape_energy,
        scaled_decay_rate_per_shell,
        injected,
    )
