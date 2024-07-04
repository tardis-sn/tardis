import logging
import numpy as np
import pandas as pd
import astropy.units as u
import radioactivedecay as rd

from tardis.energy_input.util import KEV2ERG

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def create_isotope_dicts(raw_isotope_abundance, cell_masses):
    """
    Function to create a dictionary of isotopes for each shell with their masses.

    Parameters
    ----------
    raw_isotope_abundance : pd.DataFrame
        isotope abundance in mass fractions.
    cell_masses : numpy.ndarray
        shell masses in units of g

    Returns
    -------
        isotope_dicts : Dict
            dictionary of isotopes for each shell with their ``masses``.
            Each value is abundance * cell masses.
            For eg: {0: {'Ni56': 0.1, 'Fe52': 0.2, 'Cr48': 0.3},
                    {1: {'Ni56': 0.1, 'Fe52': 0.2, 'Cr48': 0.3}} etc
    """
    isotope_dicts = {}
    for i in range(len(raw_isotope_abundance.columns)):
        isotope_dicts[i] = {}
        for (
            atomic_number,
            mass_number,
        ), abundances in raw_isotope_abundance.iterrows():
            nuclear_symbol = f"{rd.utils.Z_to_elem(atomic_number)}{mass_number}"
            isotope_dicts[i][nuclear_symbol] = (
                abundances[i] * cell_masses[i].to(u.g).value
            )

    return isotope_dicts


def create_inventories_dict(isotope_dict):
    """Function to create dictionary of inventories for each shell

    Parameters
    ----------
    isotope_dict : Dict
        dictionary of isotopes for each shell with their ``masses``.
    Returns
    -------
        inv : Dict
            dictionary of inventories for each shell
            {0: <radioactivedecay.inventory.Inventory object at 0x7f8b1c1b3d90>,
            1: <radioactivedecay.inventory.Inventory object at 0x7f8b1c1b3d90>}

    """
    inventories = {}
    for shell, isotopes in isotope_dict.items():
        inventories[shell] = rd.Inventory(isotopes, "g")

    return inventories


def calculate_total_decays(inventories, time_delta):
    """Function to create inventories of isotope for the entire simulation time.

    Parameters
    ----------
    inventories : Dict
        dictionary of inventories for each shell

    time_end : float
        End time of simulation in days.


    Returns
    -------
    cumulative_decay_df : pd.DataFrame
        total decays for x g of isotope for time 't'
    """
    time_delta = u.Quantity(time_delta, u.s)
    total_decays = {}
    for shell, inventory in inventories.items():
        total_decays[shell] = inventory.cumulative_decays(time_delta.value)

    flattened_dict = {}

    for shell, isotope_dict in total_decays.items():
        for isotope, num_of_decays in isotope_dict.items():
            new_key = isotope.replace("-", "")
            flattened_dict[(shell, new_key)] = num_of_decays

    indices = pd.MultiIndex.from_tuples(
        flattened_dict.keys(), names=["shell_number", "isotope"]
    )
    cumulative_decay_df = pd.DataFrame(
        list(flattened_dict.values()),
        index=indices,
        columns=["number_of_decays"],
    )

    return cumulative_decay_df


def create_isotope_decay_df(cumulative_decay_df, gamma_ray_lines):
    """
    Function to create a dataframe of isotopes for each shell with their decay mode, number of decays, radiation type,
    radiation energy and radiation intensity.

    Parameters
    ----------
    cumulative_decay_df : pd.DataFrame
        total decays for x g of isotope for time 't'
    gamma_ray_lines : pd.DataFrame
        gamma ray lines from nndc stored as a pandas dataframe.

    Returns
    -------
    isotope_decay_df : pd.DataFrame
        dataframe of isotopes for each shell with their decay mode, number of decays, radiation type,
        radiation energy and radiation intensity.
    """

    gamma_ray_lines = gamma_ray_lines.rename_axis(
        "isotope"
    )  # renaming "Isotope" in nndc to "isotope"
    gamma_ray_lines.drop(columns=["A", "Z"])
    gamma_ray_lines_df = gamma_ray_lines[
        ["Decay Mode", "Radiation", "Rad Energy", "Rad Intensity"]
    ]  # selecting from existing dataframe

    columns = [
        "decay_mode",
        "radiation",
        "radiation_energy_keV",
        "radiation_intensity",
    ]
    gamma_ray_lines_df.columns = columns
    isotope_decay_df = pd.merge(
        cumulative_decay_df.reset_index(),
        gamma_ray_lines_df.reset_index(),
        on=["isotope"],
    )
    isotope_decay_df = isotope_decay_df.set_index(["shell_number", "isotope"])
    isotope_decay_df["decay_mode"] = isotope_decay_df["decay_mode"].astype(
        "category"
    )
    isotope_decay_df["radiation"] = isotope_decay_df["radiation"].astype(
        "category"
    )
    isotope_decay_df["energy_per_channel_keV"] = (
        isotope_decay_df["radiation_intensity"]
        / 100.0
        * isotope_decay_df["radiation_energy_keV"]
    )
    isotope_decay_df["decay_energy_keV"] = (
        isotope_decay_df["energy_per_channel_keV"]
        * isotope_decay_df["number_of_decays"]
    )
    isotope_decay_df["decay_energy_erg"] = (
        isotope_decay_df["decay_energy_keV"] * KEV2ERG
    )

    return isotope_decay_df
