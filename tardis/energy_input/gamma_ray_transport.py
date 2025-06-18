import logging

import astropy.units as u
import numpy as np
import pandas as pd
import radioactivedecay as rd

from tardis.energy_input.energy_source import (
    get_all_isotopes,
    setup_input_energy,
)
from tardis.opacities.opacities import M_P

# Energy: keV, exported as eV for SF solver
# distance: cm
# mass: g
# time: s
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def calculate_ejecta_velocity_volume(model):
    outer_velocities = model.v_outer.to("cm/s").value
    inner_velocities = model.v_inner.to("cm/s").value
    ejecta_velocity_volume = (
        4 * np.pi / 3 * (outer_velocities**3.0 - inner_velocities**3.0)
    )

    return ejecta_velocity_volume


def calculate_average_energies(raw_isotope_abundance, gamma_ray_lines):
    """
    Function to calculate average energies of positrons and gamma rays
    from a list of gamma ray lines from nndc.

    Parameters
    ----------
    raw_isotope_abundance : pd.DataFrame
        isotope abundance in mass fractions
    gamma_ray_lines : pd.DataFrame
        decay data

    Returns
    -------
    average_energies_list : List
        list of gamma ray energies
    average_positron_energies_list : List
        list of positron energies
    gamma_ray_line_array_list : List
        list of gamma ray lines

    """
    all_isotope_names = get_all_isotopes(raw_isotope_abundance)
    all_isotope_names.sort()

    gamma_ray_line_array_list = []
    average_energies_list = []
    average_positron_energies_list = []

    gamma_ray_line_dict = {}
    average_energies = {}
    average_positron_energies = {}

    for i, isotope in enumerate(all_isotope_names):
        energy, intensity = setup_input_energy(
            gamma_ray_lines[gamma_ray_lines.index == isotope.replace("-", "")],
            "g",
        )
        average_energies_list.append(np.sum(energy * intensity))  # keV
        gamma_ray_line_array_list.append(np.stack([energy, intensity]))

        positron_energy, positron_intensity = setup_input_energy(
            gamma_ray_lines[gamma_ray_lines.index == isotope.replace("-", "")],
            "bp",
        )
        average_positron_energies_list.append(
            np.sum(positron_energy * positron_intensity)
        )

    for iso, lines in zip(all_isotope_names, gamma_ray_line_array_list):
        gamma_ray_line_dict[iso] = lines

    for iso, energy, positron_energy in zip(
        all_isotope_names, average_energies_list, average_positron_energies_list
    ):
        average_energies[iso] = energy
        average_positron_energies[iso] = positron_energy

    return (
        average_energies,
        average_positron_energies,
        gamma_ray_line_dict,
    )


def calculate_average_power_per_mass(energy_per_mass, time_delta):
    # Time averaged energy per mass for constant packet count
    average_power_per_mass = energy_per_mass / (time_delta)

    return average_power_per_mass


def iron_group_fraction_per_shell(model):
    # Taking iron group to be elements 21-30
    # Used as part of the approximations for photoabsorption and pair creation
    # Dependent on atomic data
    return model.abundance.loc[(21):(30)].sum(axis=0)
