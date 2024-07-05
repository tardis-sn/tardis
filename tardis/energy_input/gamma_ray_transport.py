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


def get_nuclide_atomic_number(nuclide):
    return rd.Nuclide(nuclide).Z


def get_chain_decay_power_per_ejectamass(
    inventory,
    time,
    time_start,
    average_energies,
    average_positron_energies,
    tau_end,
):
    # total decay power per mass [erg/s/g] for a given decaypath
    # only decays at the end of the chain contributed from the initial abundance of the top of the chain are counted
    # (these can be can be same for a chain of length one)

    # need to find a way to connect inventory_start with inventory_end to obtain start and end of chains
    # hopefully can just do it all at once

    # start_isotope is e.g. Ni56 while end_isotope would be e.g. Fe56
    inventory_end = inventory.decay(time - time_start)
    # contribution to the end nuclide abundance from the top of chain (could be a length-one chain Z,A_top = Z,A_end
    # so contribution would be from init abundance only)

    decaypower = []
    progeny = inventory.progeny()

    for start_isotope in progeny:
        # ignore branching probs for now
        end_isotope = progeny[start_isotope][0]

        if rd.Nuclide(end_isotope).half_life("readable") == "stable":
            end_isotope = start_isotope

        # i.e. need to get end_isotope(s) by finding out which ones have no progeny left
        # this should be the total number of "end chain" isotopes after (time - time_start)
        endnucabund = inventory_end.mass_fractions()[end_isotope]

        # this is the total decay energy from gamma-rays and positrons for the end of chain isotope
        # endecay = get_decaypath_lastnucdecayenergy(decaypathindex)
        endecay = (
            average_energies[end_isotope]
            + average_positron_energies[end_isotope]
        )

        print("Decay energy, abundance, tau")
        print(endecay)
        print(endnucabund)
        print(tau_end)

        decaypower.append(
            endecay
            * endnucabund
            / tau_end
            / (rd.Nuclide(start_isotope).atomic_mass * M_P)
        )

    return decaypower


def calculate_ejecta_velocity_volume(model):
    outer_velocities = model.v_outer.to("cm/s").value
    inner_velocities = model.v_inner.to("cm/s").value
    ejecta_velocity_volume = (
        4 * np.pi / 3 * (outer_velocities**3.0 - inner_velocities**3.0)
    )

    return ejecta_velocity_volume


def calculate_total_decays_old(inventories, time_delta):
    """Function to create inventories of isotope

    Parameters
    ----------
    model : tardis.Radial1DModel
        The tardis model to calculate gamma ray propagation through

    time_end : float
        End time of simulation in days

    Returns
    -------
        Total decay list : List
            list of total decays for x g of isotope for time 't'
    """
    time_delta = u.Quantity(time_delta, u.s)
    total_decays = {}
    for shell, isotopes in inventories.items():
        total_decays[shell] = {}
        for isotope, name in isotopes.items():
            # decays = name.decay(time_delta.value, "s")
            total_decays[shell][isotope] = name.cumulative_decays(
                time_delta.value
            )
    return total_decays


def create_isotope_dicts_old(raw_isotope_abundance, cell_masses):
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
            isotope_dicts[i][atomic_number, mass_number] = {}
            nuclear_symbol = f"{rd.utils.Z_to_elem(atomic_number)}{mass_number}"
            isotope_dicts[i][atomic_number, mass_number][nuclear_symbol] = (
                abundances[i] * cell_masses[i].to(u.g).value
            )

    return isotope_dicts


def create_inventories_dict_old(isotope_dict):
    """Function to create dictionary of inventories for each shell

    Parameters
    ----------
    isotope_dict : Dict
        dictionary of isotopes for each shell with their ``masses``.

    Returns
    -------
        inv : Dict
            dictionary of inventories for each shell
            For eg: {0: {'Ni56': <radioactivedecay.Inventory at 0x7f7d2c0d0e50>,
                         'Co56': <radioactivedecay.Inventory at 0x7f7d2c0d0e50>},
                    {1: {'Ni56': <radioactivedecay.Inventory at 0x7f7d2c0d0e50>,
                         'Co56': <radioactivedecay.Inventory at 0x7f7d2c0d0e50>}} etc
    """
    inv = {}
    for shell, isotopes in isotope_dict.items():
        inv[shell] = {}
        for isotope, name in isotopes.items():
            inv[shell][isotope] = rd.Inventory(name, "g")

    return inv


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


def get_taus(raw_isotope_abundance):
    """
    Function to calculate taus for each isotope

    Parameters
    ----------
    raw_isotope_abundance : pd.DataFrame
        isotope abundance in mass fractions

    Returns
    -------
    taus : Dict
        dictionary of taus for each isotope
    parents : Dict
        dictionary of parents for each isotope
    """
    all_isotope_names = get_all_isotopes(raw_isotope_abundance)
    all_isotope_names.sort()

    taus = {}
    parents = {}
    for isotope in all_isotope_names:
        nuclide = rd.Nuclide(isotope)
        taus[isotope] = nuclide.half_life() / np.log(2)
        child = nuclide.progeny()
        if child is not None:
            for c in child:
                if rd.Nuclide(c).half_life("readable") != "stable":
                    # this is a dict of child: parent intended to find
                    # the parents of a given isotope.
                    # if there is no parent, there is no item.
                    parents[c] = isotope

    return taus, parents


def decay_chain_energies(
    average_energies,
    total_decays,
):
    """
    Function to calculate decay chain energies.

    Parameters
    ----------
    raw_isotope_abundance : pd.DataFrame
        isotope abundance in mass fractions
    average_energies_list : List
        list of gamma ray energies
    average_positron_energies_list : List
        list of positron energies
    gamma_ray_line_array_list : List
        list of gamma ray lines
    total_decays : Dict
        dictionary of total decays for each isotope in each shell

    Returns
    -------
    decay_energy : Dict
        dictionary of decay chain energies for each isotope in each shell

    """
    decay_energy = {}
    for shell, isotopes in total_decays.items():
        decay_energy[shell] = {}
        for name, isotope in isotopes.items():
            decay_energy[shell][name] = {}
            for iso, dps in isotope.items():
                decay_energy[shell][name][iso] = dps * average_energies[iso]

    return decay_energy


def fractional_decay_energy(decay_energy):
    """Function to calculate fractional decay energy

    Parameters
    ----------
    decay_energy : Dict
        dictionary of decay chain energies for each isotope in each shell

    Returns
    -------
    fractional_decay_energy : Dict
        dictionary of fractional decay chain energies for each isotope in each shell
    """
    fractional_decay_energy = {
        shell: {
            parent_isotope: {
                isotopes: (
                    decay_energy[shell][parent_isotope][isotopes]
                    / sum(decay_energy[shell][parent_isotope].values())
                    if decay_energy[shell][parent_isotope][isotopes] != 0.0
                    else 0.0
                )
                for isotopes in decay_energy[shell][parent_isotope]
            }
            for parent_isotope in decay_energy[shell]
        }
        for shell in decay_energy
    }

    return fractional_decay_energy


def calculate_energy_per_mass(decay_energy, raw_isotope_abundance, cell_masses):
    """
    Function to calculate decay energy per mass for each isotope chain.

    Parameters
    ----------
    decay_energy : Dict
        dictionary of decay chain energies for each isotope in each shell
    raw_isotope_abundance : pd.DataFrame
        isotope abundance in mass fractions.
    cell_masses : numpy.ndarray
        shell masses in units of g

    Returns
    -------
    energy_per_mass : pd.DataFrame
        decay energy per mass for each isotope chain
        For e.g Ni56 has 2 decay chains:
        Ni56 -> Co56 -> Fe56. It will calculate the decay energy per mass for each chain
        and store as a dataframe.
    """
    energy_dict = {}
    for shell, isotopes in decay_energy.items():
        energy_dict[shell] = {}
        for name, isotope in isotopes.items():
            energy_dict[shell][name] = sum(isotope.values())

    energy_list = []
    for shell, isotopes in energy_dict.items():
        for isotope, energy in isotopes.items():
            energy_list.append(
                {
                    "shell": shell,
                    "atomic_number": isotope[0],
                    "mass_number": isotope[1],
                    "value": energy,
                }
            )

    df = pd.DataFrame(energy_list)
    energy_df = pd.pivot_table(
        df,
        values="value",
        index=["atomic_number", "mass_number"],
        columns="shell",
    )

    energy_per_mass = energy_df.divide(
        (raw_isotope_abundance * cell_masses).to_numpy(), axis=0
    )

    return energy_per_mass, energy_df


def distribute_packets(decay_energy, total_energy, num_packets):
    packets_per_isotope = {}
    for shell, isotopes in decay_energy.items():
        packets_per_isotope[shell] = {}
        for name, isotope in isotopes.items():
            packets_per_isotope[shell][name] = {}
            for line, energy in isotope.items():
                packets_per_isotope[shell][name][line] = int(
                    energy / total_energy * num_packets
                )

    packets_per_isotope_list = []
    for shell, parent_isotope in packets_per_isotope.items():
        for isotopes, isotope_dict in parent_isotope.items():
            for name, value in isotope_dict.items():
                packets_per_isotope_list.append(
                    {
                        "shell": shell,
                        "element": name,
                        "value": value,
                    }
                )

    df = pd.DataFrame(packets_per_isotope_list)
    packets_per_isotope_df = pd.pivot_table(
        df,
        values="value",
        index="element",
        columns="shell",
    )

    return packets_per_isotope_df


def packets_per_isotope(fractional_decay_energy, decayed_packet_count_dict):
    packets_per_isotope = {
        shell: {
            parent_isotope: {
                isotopes: fractional_decay_energy[shell][parent_isotope][
                    isotopes
                ]
                * decayed_packet_count_dict[shell][parent_isotope]
                for isotopes in fractional_decay_energy[shell][parent_isotope]
            }
            for parent_isotope in fractional_decay_energy[shell]
        }
        for shell in fractional_decay_energy
    }

    packets_per_isotope_list = []
    for shell, parent_isotope in packets_per_isotope.items():
        for isotopes, isotope_dict in parent_isotope.items():
            for name, value in isotope_dict.items():
                packets_per_isotope_list.append(
                    {
                        "shell": shell,
                        "element": name,
                        "value": value,
                    }
                )

    df = pd.DataFrame(packets_per_isotope_list)
    packets_per_isotope_df = pd.pivot_table(
        df,
        values="value",
        index="element",
        columns="shell",
    )

    return packets_per_isotope_df


def calculate_average_power_per_mass(energy_per_mass, time_delta):
    # Time averaged energy per mass for constant packet count
    average_power_per_mass = energy_per_mass / (time_delta)

    return average_power_per_mass


def iron_group_fraction_per_shell(model):
    # Taking iron group to be elements 21-30
    # Used as part of the approximations for photoabsorption and pair creation
    # Dependent on atomic data
    return model.abundance.loc[(21):(30)].sum(axis=0)
