import logging
import numpy as np
import pandas as pd
import astropy.units as u
import radioactivedecay as rd


from tardis.energy_input.energy_source import (
    get_nuclear_lines_database,
)

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
            For eg: {0: {(28, 56): {'Ni56': 0.0001}, (27, 57): {'Co56': 0.0001}}
                    {1: {(28, 56): {'Ni56': 0.0001}, (27, 57): {'Co56': 0.0001}}} etc

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
            For eg: {0: {'Ni56': <radioactivedecay.Inventory>,
                         'Co56': <radioactivedecay.Inventory>},
                    {1: {'Ni56': <radioactivedecay.Inventory>,
                         'Co56': <radioactivedecay.Inventory>}} etc
    """
    inv = {}
    for shell, isotopes in isotope_dict.items():
        inv[shell] = rd.Inventory(isotopes, "g")

    return inv


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
    for shell, isotopes in inventories.items():
        total_decays[shell] = isotopes.cumulative_decays(time_delta.value)

    flattened_dict = {}

    for shell, isotope_dict in total_decays.items():
        for isotope, decay_value in isotope_dict.items():
            new_key = isotope.replace("-", "")
            flattened_dict[(shell, new_key)] = decay_value

    indices = pd.MultiIndex.from_tuples(
        flattened_dict.keys(), names=["Shell_number", "Isotope"]
    )
    cumulative_decay_df = pd.DataFrame(
        list(flattened_dict.values()),
        index=indices,
        columns=["Number_of_decays"],
    )

    return cumulative_decay_df


def create_isotope_decay_df(cumulative_decay_df, atom_data_file):
    """
    Function to create a dataframe of isotopes for each shell with their decay mode, number of decays, radiation type,
    radiation energy and radiation intensity.

    Parameters
    ----------
    cumulative_decay_df : pd.DataFrame
        dataframe of isotopes for each shell with their number of decays.
    atom_data_file : str
        path to the atomic data file.

    Returns
    -------
    isotope_decay_df : pd.DataFrame
        dataframe of isotopes for each shell with their decay mode, number of decays, radiation type,
        radiation energy and radiation intensity.
    """

    gamma_ray_lines = get_nuclear_lines_database(atom_data_file)
    gamma_ray_lines.reset_index(inplace=True)
    gamma_ray_lines.set_index(["Isotope", "A", "Z"], inplace=True)

    gamma_ray_lines_df = gamma_ray_lines[
        ["Decay Mode", "Radiation", "Rad Energy", "Rad Intensity"]
    ]

    isotope_decay_df = pd.merge(
        cumulative_decay_df.reset_index(),
        gamma_ray_lines_df.reset_index(),
        on=["Isotope"],
    )
    isotope_decay_df
    isotope_decay_df.set_index(["Shell_number", "Isotope"], inplace=True)

    return isotope_decay_df
