import astropy.units as u
import tardis.constants as const
import numpy as np
import scipy.linalg as linalg

M_E = const.m_e.to_cgs().value
H_BAR = const.hbar.to_cgs().value
E_CHARGE = const.e.to_cgs().value

# Eqn 7 set up energy grid, bin-wise integration, multiply S(E) by the energy grid


def coulomb_loss_function(energy, electron_number_density, number_density):
    """Calculates the Coulomb loss function for a given energy,
    electron number density, and atom number density

    Parameters
    ----------
    energy : float
        electron energy
    electron_number_density : float

    number_density : float

    Returns
    -------
    float
        Coulomb loss energy
    """
    plasma_frequency = 56414.6 * np.sqrt(electron_number_density)
    zeta_electron = 2.0 * H_BAR * plasma_frequency
    electron_fraction = electron_number_density / number_density

    if energy > 14:
        assert 2 * energy > zeta_electron
        return (
            electron_fraction
            * (2 * np.pi * E_CHARGE ** 4)
            / energy
            * np.log(2 * energy / zeta_electron)
        )
    else:
        v = np.sqrt(2 * energy / M_E)
        euler_gamma = 0.577215664901532
        return (
            electron_fraction
            * (2 * np.pi * E_CHARGE ** 4)
            / energy
            * np.log(
                M_E * v ** 3 / (euler_gamma * E_CHARGE ** 2 * plasma_frequency)
            )
        )


def cross_section(
    energy, initial_electron_energy, oscillator_strength, van_regemorter_fit
):
    # Probably more useful from macro atom
    hydrogen_ionization_potential = 1

    k = initial_electron_energy / 13.60

    return (
        (8 * np.pi)
        / np.sqrt(3)
        * (1 / k ** 2)
        * (hydrogen_ionization_potential / energy)
        * oscillator_strength
        * van_regemorter_fit
        * const.a0 ** 2
    )


def collisional_cross_section(energy):
    # Younger 1981 apparently
    # But also used in macro atom
    return False


def ar_xs(energy, ion_potential, A, B, C, D):
    u = energy / ion_potential
    if u <= 1:
        return 0

    return (
        1e-14
        * (
            A * (1 - 1 / u)
            + B * (1 - 1 / u) ** 2
            + C * np.log(u)
            + D * np.log(u) / u
        )
        / (u * ion_potential ** 2)
    )


def get_arxs_array_shell(energy_grid, shell):
    ar_xs_array = np.array(
        [
            ar_xs(
                energy, shell.ion_potential, shell.A, shell.B, shell.C, shell.D
            )
            for energy in energy_grid
        ]
    )

    return ar_xs_array


def get_J(Z, ionization_stage, ion_potential):
    # returns an energy in eV
    # values from Opal et al. 1971 as applied by Kozma & Fransson 1992
    if ionization_stage == 1:
        if Z == 2:  # He I
            return 15.8
        elif Z == 10:  # Ne I
            return 24.2
        elif Z == 18:  # Ar I
            return 10.0

    return 0.6 * ion_potential


def get_index(indexed_energy, energy_grid):
    assert indexed_energy >= energy_grid[0]
    assert indexed_energy < (
        energy_grid[-1] + (energy_grid[1] - energy_grid[0])
    )

    for i, energy in enumerate(energy_grid):
        if energy < indexed_energy:
            index = i

    return index


def spencer_fano_matrix_add_ionization_shell(
    energy_grid, delta_energy, points, number_ion, shell, spencer_fano_matrix
):
    """contains the terms related to ionisation cross sections"""
    ion_potential = shell.ion_potential
    J = get_J(shell.Z, shell.ionization_stage, ion_potential)

    ar_xs_array = get_arxs_array_shell(energy_grid, shell)

    if ion_potential <= energy_grid[0]:
        xs_start_index = 0
    else:
        xs_start_index = get_index(ion_potential, energy_grid)

    for i, energy in enumerate(energy_grid):
        # // endash ranges from en to SF_EMAX, but skip over the zero-cross section points
        j_start = i if i > xs_start_index else xs_start_index
        if 2 * energy + ion_potential < energy_grid[-1] + (
            energy_grid[1] - energy_grid[0]
        ):
            secondintegralstartindex = get_index(
                2 * energy + ion_potential, energy_grid
            )
        else:
            secondintegralstartindex = points + 1

        # integral/J of 1/[1 + (epsilon - ion_potential) / J] for epsilon = en + ion_potential
        for j in range(j_start, points):
            # j is the matrix column index which corresponds to the piece of the
            # integral at y(E') where E' >= E and E' = envec(j)
            energy_dash = energy_grid[j]
            prefactor = (
                nnion
                * ar_xs_array[j]
                / np.atan((energy_dash - ion_potential) / 2.0 / J)
                * delta_energy
            )
            assert not np.isnan(prefactor)
            assert not np.isinf(prefactor)
            # assert prefactor >= 0

            # J * atan[(epsilon - ionpot_ev) / J] is the indefinite integral of
            # 1/(1 + (epsilon - ionpot_ev)^2/ J^2) d_epsilon
            # in Kozma & Fransson 1992 equation 4

            # KF 92 limit
            epsilon_upper = (energy_dash + ion_potential) / 2
            # Li+2012 limit
            # epsilon_upper = (endash + en) / 2

            int_eps_upper = np.arctan((epsilon_upper - ion_potential) / J)

            epsilon_lower = energy_dash - energy
            int_eps_lower = np.arctan((epsilon_lower - ion_potential) / J)

            spencer_fano_matrix[i, j] += prefactor * (
                int_eps_upper - int_eps_lower
            )

            epsilon_lower = energy + ion_potential
            epsilon_upper = (energy_dash + ion_potential) / 2
            # endash ranges from 2 * en + ionpot_ev to SF_EMAX
            if j >= secondintegralstartindex + 1:
                # int_eps_upper = atan((epsilon_upper - ionpot_ev) / J)
                int_eps_lower = np.arctan((epsilon_lower - ion_potential) / J)
                if epsilon_lower > epsilon_upper:
                    print(
                        j,
                        secondintegralstartindex,
                        epsilon_lower,
                        epsilon_upper,
                    )
                assert epsilon_lower <= epsilon_upper

                spencer_fano_matrix[i, j] -= prefactor * (
                    int_eps_upper - int_eps_lower
                )

    return spencer_fano_matrix


def solve_spencer_fano(
    energy_grid,
    source_vector,
    electron_number_density,
    number_density,
    ions,
    deposition_density,
    ion_coll_data,
):
    delta_energy = energy_grid[1] - energy_grid[0]
    points = len(energy_grid)

    initial_energy = np.dot(energy_grid, source_vector) * delta_energy

    # set up the constant vector
    # sum of source vector times dE
    # same size as energy grid
    constant_vector = np.zeroes(points)
    for i in range(points):
        for j in range(i, points):
            constant_vector[i] += source_vector[j] * delta_energy

    spencer_fano_matrix = np.zeroes((points, points))
    for i in range(points):
        energy = energy_grid[i]
        spencer_fano_matrix[i, i] += coulomb_loss_function(
            energy, electron_number_density, number_density
        )

    df_transitions = {}

    for Z, ion_stage in ions:
        number_ion = ion_population(Z, ion_stage)

        for index, shell in ion_coll_data:
            assert shell.ion_potential >= energy_grid[0]
            spencer_fano_matrix = spencer_fano_matrix_add_ionization_shell(
                energy_grid,
                delta_energy,
                points,
                number_ion,
                shell,
                spencer_fano_matrix,
            )

        if not no_excitation:
            # hellish madness
            pass

    lu_and_piv = linalg.lu_factor(spencer_fano_matrix, overwrite_a=False)
    y_vector_reference = linalg.lu_solve(lu_and_piv, constant_vector, trans=0)
    y_vector = y_vector_reference * deposition_density / initial_energy

    return y_vector, df_transitions


def electron_spectrum(energy):
    """number of electrons at energy E

    Parameters
    ----------
    energy : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    return False


def heating_fraction(
    energy,
    mean_initial_energy,
    electron_number_density,
    number_density,
    lowest_energy,
    max_energy,
):
    """Calculates the heating fraction of electrons

    Parameters
    ----------
    energy : float

    mean_initial_energy : float

    electron_number_density : float

    number_density : float

    lowest_energy : float

    max_energy : float

    Returns
    -------
    float
        fraction of energy going into heating
    """
    mean_fraction = 1 / mean_initial_energy

    constant = (
        lowest_energy
        * electron_spectrum(lowest_energy)
        * coulomb_loss_function(
            lowest_energy, electron_number_density, number_density
        )
    )

    degradation_energies = np.linspace(lowest_energy, max_energy)
    degradation_spectrum = electron_spectrum(degradation_energies)
    coulomb_loss_spectrum = coulomb_loss_function(
        degradation_energies, electron_number_density, number_density
    )

    integral_degradation = np.trapz(
        degradation_spectrum * coulomb_loss_spectrum
    )

    number_of_things = 1 * 1

    integral_number = np.trapz(number_of_things)

    return (
        mean_fraction * integral_degradation
        + mean_fraction * constant
        + mean_fraction * integral_number
    )


def excitation_fraction(
    energy, max_energy, number_density, mean_initial_energy
):
    """Excitation fraction of electron energy

    Parameters
    ----------
    energy : float
        [description]
    max_energy : float
        [description]
    number_density : float
        [description]
    mean_initial_energy : float
        [description]

    Returns
    -------
    float
        fraction of energy going into excitation
    """
    degradation_energies = np.linspace(energy, max_energy)
    degradation_spectrum = electron_spectrum(degradation_energies)

    cross_sections = cross_section(degradation_energies)

    integral_degradation = np.trapz(degradation_spectrum * cross_sections)

    return (
        (number_density * energy) / mean_initial_energy * integral_degradation
    )


def ionization_fraction(
    energy, max_energy, number_density, mean_initial_energy
):
    """Ionization fraction of electron energy

    Parameters
    ----------
    energy : float
        [description]
    max_energy : float
        [description]
    number_density : float
        [description]
    mean_initial_energy : float
        [description]

    Returns
    -------
    float
        Ionization fraction of electron energy
    """
    degradation_energies = np.linspace(energy, max_energy)
    degradation_spectrum = electron_spectrum(degradation_energies)

    q = collisional_cross_section(degradation_energies)

    integral_degradation = np.trapz(degradation_spectrum * q)

    return (
        (number_density * energy) / mean_initial_energy * integral_degradation
    )
