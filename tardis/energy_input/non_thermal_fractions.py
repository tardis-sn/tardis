import astropy.units as u
import tardis.constants as const
import numpy as np
import scipy.linalg as linalg

import spencer_fano as sf

M_E = const.m_e.to_cgs().value
H_BAR = const.hbar.to_cgs().value
E_CHARGE = const.e.to_cgs().value


def probability_secondary_electron(
    primary_electron_energy,
    ionization_potential,
    J,
    secondary_electron_energy=-1,
    epsilon=-1,
):
    assert secondary_electron_energy >= 0 or epsilon >= 0
    if secondary_electron_energy < 0:
        secondary_electron_energy = epsilon - ionization_potential
    if epsilon < 0:
        epsilon = secondary_electron_energy + ionization_potential
    return (
        1.0
        / J
        / np.arctan((primary_electron_energy - ionization_potential) / 2.0 / J)
        / (1 + ((secondary_electron_energy / J) ** 2))
    )


def calculate_number_electron(energy, energy_grid, electron_spectrum):
    # Kozma & Fransson equation 6.
    # Something related to a number of electrons, needed to calculate the heating fraction in equation 3
    # not valid for energy > E_0
    if energy == 0.0:
        return 0.0

    # Andreas cached values

    N_e = 0.0

    for atomic_number, ion_number in ions:
        N_e_ion = 0.0
        nnion = ionpopdict[(atomic_number, ion_number)]
        dfcollion_thision = dfcollion.query(
            "Z == @Z and ionstage == @ionstage", inplace=False
        )

        for index, shell in dfcollion_thision.iterrows():
            ionization_potential = shell.ionpot_ev

            energy_lambda = min(
                energy_grid.energy_max - energy, energy + ionization_potential
            )
            J = sf.get_J(shell.Z, shell.ionstage, ionization_potential)

            arnaud_cross_section_array = (
                sf.get_arnaud_cross_section_array_shell(energy_grid.grid, shell)
            )

            # integral from ionpot to energy_lambda

            delta_energy_dash = (energy_lambda - ionization_potential) / 1000.0
            if delta_energy_dash >= 0:
                endashlist = np.arange(
                    ionization_potential, energy_lambda, delta_energy_dash
                )
                for energy_dash in endashlist:
                    i = sf.get_index(energy + energy_dash, energy_grid.grid)
                    N_e_ion += (
                        electron_spectrum[i]
                        * arnaud_cross_section_array[i]
                        * probability_secondary_electron(
                            primary_electron_energy=energy + energy_dash,
                            epsilon=energy_dash,
                            ionization_potential=ionization_potential,
                            J=J,
                        )
                        * delta_energy_dash
                    )

            # // integral from 2E + I up to E_max
            delta_energy_dash = (
                energy_grid.energy_max - (2 * energy + ionization_potential)
            ) / 100.0
            if delta_energy_dash >= 0:
                endashlist = np.arange(
                    2 * energy + ionization_potential,
                    energy_grid.energy_max,
                    delta_energy_dash,
                )
                for energy_dash in endashlist:
                    i = sf.get_index(energy_dash, energy_grid.grid)
                    N_e_ion += (
                        electron_spectrum[i]
                        * arnaud_cross_section_array[i]
                        * probability_secondary_electron(
                            primary_electron_energy=energy_dash,
                            epsilon=energy_dash,
                            ionization_potential=ionization_potential,
                            J=J,
                        )
                        * delta_energy_dash
                    )

        N_e += nnion * N_e_ion

    for Z, ion_stage in ions:
        for _, row in dftransitions[(Z, ion_stage)].iterrows():
            nnlevel = row.lower_pop
            epsilon_trans_ev = row.epsilon_trans_ev
            if epsilon_trans_ev >= energy_grid.energy_min:
                i = sf.get_index(energy + epsilon_trans_ev, energy_grid.grid)
                cross_section_excitation_vector = (
                    sf.excitation_cross_section_vector(energy_grid, row)
                )
                N_e += (
                    nnlevel
                    * epsilon_trans_ev
                    * cross_section_excitation_vector[i]
                    * electron_spectrum[i]
                )

    # source term not here because it should be zero at the low end anyway

    return N_e


def heating_fraction(
    energy_grid,
    electron_spectrum,
    electron_number_density,
    number_density,
    energy_deposition_density,
):
    """Calculates the heating fraction of electrons

    Parameters
    ----------
    Returns
    -------
    float
        fraction of energy going into heating
    """

    for i, energy in enumerate(energy_grid.grid):
        weight = 1 if (i == 0 or i == energy_grid.size - 1) else 2
        fraction_heating += (
            0.5
            * weight
            * sf.coulomb_loss_function(
                energy, electron_number_density, number_density
            )
            * electron_spectrum[i]
            * energy_grid.delta_energy
            / energy_deposition_density,
        )

    fraction_heating += (
        energy_grid.energy_min
        * electron_spectrum[0]
        * sf.coulomb_loss_function(
            energy_grid.energy_min, electron_number_density, number_density
        )
        / energy_deposition_density
    )

    fraction_heating_number_electron = 0.0
    delta_energy = energy_grid.energy_min / 10.0  # seems arbitrary
    for energy in np.arange(0.0, energy_grid.energy_min, delta_energy):
        number_electron = calculate_number_electron(
            energy,
            energy_grid,
            ions,
            ion_population,
            collision_data,
            electron_spectrum,
            transitions_dict,
        )

        fraction_heating_number_electron += (
            number_electron * energy * delta_energy / energy_deposition_density
        )

    fraction_heating += fraction_heating_number_electron
    return fraction_heating


def excitation_fraction_per_ion(
    energy_grid, transitions_dict, electron_spectrum, energy_deposition_density
):
    """Excitation fraction of electron energy

    Parameters
    ----------

    Returns
    -------
    float
        fraction of energy going into excitation
    """
    cross_section_excitation_vector_sum = np.zeros(energy_grid.size)

    for _, row in transitions_dict.iterrows():
        number_in_level = row.lower_pop
        cross_section_excitation_vector_sum += (
            number_in_level
            * row.epsilon_trans_ev
            * sf.excitation_cross_section_vector(energy_grid, row)
        )

    return (
        np.dot(cross_section_excitation_vector_sum, electron_spectrum)
        * energy_grid.delta_energy
        / energy_deposition_density
    )


def ionization_fraction_per_ion(
    energy_grid,
    electron_spectrum,
    energy_deposition_density,
    ion_collision_data,
    ion_number_density,
):
    """Ionization fraction of electron energy per ion

    Parameters
    ----------

    Returns
    -------
    float
        Ionization fraction of electron energy
    """
    fractional_ionization_per_ion = 0.0
    for index, shell in ion_collision_data.iterrows():
        cross_section_start_index = sf.get_index(
            shell.ion_potential, energy_grid.grid
        )
        arnaud_cross_section_array = sf.get_arnaud_cross_section_array_shell(
            energy_grid.grid, shell
        )

        fractional_ionization_of_shell = (
            ion_number_density
            * shell.ionpot_ev
            * np.dot(electron_spectrum, arnaud_cross_section_array)
            * energy_grid.delta_energy
            / energy_deposition_density
        )

        if fractional_ionization_of_shell > 1:
            fractional_ionization_of_shell = 0.0
            print(
                "Ignoring fractional_ionization_of_shell of {}.".format(
                    fractional_ionization_of_shell
                )
            )
        fractional_ionization_per_ion += fractional_ionization_of_shell
        print(
            "  cross section at {:2.7f} eV and {} eV {:.2e} and {:.2e}".format(
                energy_grid.grid[cross_section_start_index + 1],
                energy_grid.energy_max,
                arnaud_cross_section_array[cross_section_start_index + 1],
                arnaud_cross_section_array[-1],
            )
        )

    return fractional_ionization_per_ion


def total_fractions(
    energy_grid,
    electron_spectrum,
    transitions_dict,
    ion_collision_data,
    ion_populations,
    energy_deposition_density,
):
    fraction_ionization = 0.0
    fraction_excitation = 0.0

    for atomic_number, ion_number in ion_populations:
        ion_number_density = ion_populations[(atomic_number, ion_number)]

        ion_collision_data_current = ion_collision_data.query("This ion")

        fraction_ionization += ionization_fraction_per_ion(
            energy_grid,
            electron_spectrum,
            energy_deposition_density,
            ion_collision_data_current,
            ion_number_density,
        )

        fraction_excitation += excitation_fraction_per_ion(
            energy_grid,
            transitions_dict[(atomic_number, ion_number)],
            electron_spectrum,
            energy_deposition_density,
        )

    fraction_heating = heating_fraction(electron_number_density, number_density)

    return fraction_ionization, fraction_excitation, fraction_heating
