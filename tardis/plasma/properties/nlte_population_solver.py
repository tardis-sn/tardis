import pandas as pd
import numpy as np
import logging
from scipy.linalg import lu_factor, lu_solve

from tardis.plasma.properties.base import ProcessingPlasmaProperty

__all__ = [
    "NLTEPopulationSolver",
]

logger = logging.getLogger(__name__)

NLTE_POPULATION_SOLVER_MAX_ITERATIONS = 100
NLTE_POPULATION_SOLVER_TOLERANCE = 1e-3


class NLTEPopulationSolver(ProcessingPlasmaProperty):
    outputs = ("ion_number_density", "electron_densities")

    def calculate(
        self,
        gamma,
        alpha_sp,
        alpha_stim,
        coll_ion_coeff,
        coll_recomb_coeff,
        partition_function,
        levels,
        level_boltzmann_factor,
        phi,
        rate_matrix_index,
        number_density,
        nlte_excitation_species,
    ):
        """Calculates ion number densities and electron densities using NLTE ionization using an iterative LU decomposition approach.

        Parameters
        ----------
        gamma : pandas.DataFrame
            The rate coefficient for radiative ionization.
        alpha_sp : pandas.DataFrame
            The rate coefficient for spontaneous recombination.
        alpha_stim : pandas.DataFrame
            The rate coefficient for stimulated recombination.
        coll_ion_coeff : pandas.DataFrame
            The rate coefficient for collisional ionization in the Seaton
            approximation.
        coll_recomb_coeff : pandas.DataFrame
            The rate coefficient for collisional recombination.
        partition_function : pandas.DataFrame
            General partition function. Indexed by atomic number, ion number.
        levels : MultiIndex
            (atomic_number, ion_number, level_number)
            Index of filtered atomic data.
        level_boltzmann_factor : pandas.DataFrame
            General Boltzmann factor.
        phi : pandas.DataFrame
            Saha Factors.
        rate_matrix_index : MultiIndex
            (atomic_number, ion_number, treatment type)
            If ion is treated in LTE or nebular ionization, 3rd index is "lte_ion",
            if treated in NLTE ionization, 3rd index is "nlte_ion".
        number_density : pandas.DataFrame
            Number density in each shell for each species.
        nlte_excitation_species : list
            List of species treated in NLTE excitation.

        Returns
        -------
        ion_number_densities : pandas.DataFrame
            Number density with NLTE ionization treatment.
        electron_densities : Series
            Electron density with NLTE ionization treatment.
        """

        (
            total_photo_ion_coefficients,
            total_rad_recomb_coefficients,
            total_coll_ion_coefficients,
            total_coll_recomb_coefficients,
        ) = self.prepare_ion_recomb_coefficients_nlte_ion(
            gamma,
            alpha_sp,
            alpha_stim,
            coll_ion_coeff,
            coll_recomb_coeff,
            partition_function,
            levels,
            level_boltzmann_factor,
        )

        # TODO: Don't create the rate_matrix_index with n_e in the first place
        rate_matrix_index = rate_matrix_index.drop("n_e")

        initial_electron_densities = number_density.sum(axis=0)
        atomic_numbers = rate_matrix_index.get_level_values(
            "atomic_number"
        ).unique()

        index = rate_matrix_index.droplevel("level_number")
        electron_densities = initial_electron_densities.copy()
        ion_number_density = pd.DataFrame(0.0, index=index, columns=phi.columns)
        ion_numbers = index.get_level_values("ion_number")

        # Ordering of loops is important to allow for
        # easier parallelization in the future
        logger.info("Starting NLTE ionization solver")
        for shell in phi.columns:
            iteration = 0
            balance_vector = self.calculate_balance_vector(
                number_density[shell],
                rate_matrix_index,
            )
            while iteration < NLTE_POPULATION_SOLVER_MAX_ITERATIONS:
                rate_matrix = self.calculate_rate_matrix(
                    atomic_numbers,
                    phi[shell],
                    electron_densities[shell],
                    rate_matrix_index,
                    total_photo_ion_coefficients[shell],
                    total_rad_recomb_coefficients[shell],
                    total_coll_ion_coefficients[shell],
                    total_coll_recomb_coefficients[shell],
                )
                # TODO: Solve for each element individually
                # and handle errors in the solver
                ion_solution = self.solve(rate_matrix, balance_vector)
                electron_solution = np.sum(ion_numbers * ion_solution)

                assert np.all(
                    ion_solution >= 0.0
                ), "Negative ion number density found, solver failed."
                assert (
                    electron_solution >= 0.0
                ), "Negative electron density found, solver failed."

                delta_ion, delta_electron = self.calculate_delta(
                    ion_solution,
                    electron_solution,
                    ion_number_density[shell],
                    electron_densities[shell],
                )
                ion_number_density[shell] = ion_solution
                electron_densities[shell] = electron_solution
                if (
                    np.all(np.abs(delta_ion) < NLTE_POPULATION_SOLVER_TOLERANCE)
                    and np.abs(delta_electron)
                    < NLTE_POPULATION_SOLVER_TOLERANCE
                ):
                    ion_number_density[shell] = ion_solution
                    electron_densities[shell] = electron_solution
                    logger.debug(
                        "NLTE ionization solver converged after {} iterations for shell {}".format(
                            iteration, shell
                        )
                    )
                    break

                iteration += 1
                if iteration == NLTE_POPULATION_SOLVER_MAX_ITERATIONS:
                    logger.warning(
                        "NLTE ionization solver did not converge for shell {} ".format(
                            shell
                        )
                    )
                    break

        logger.info("NLTE ionization solver finished")

        return ion_number_density, electron_densities

    @staticmethod
    def calculate_delta(
        ion_solution, electron_solution, ion_number_density, electron_densities
    ):
        """Calculates relative change between new solution and old value.

        Parameters
        ----------
        ion_solution : numpy.array
            Solution vector for the NLTE ionization solver.
        electron_solution : float
            Solution for the electron density.

        Returns
        -------
        numpy.array
            Change in ion number densities.
        float
            Change in electron density.
        """
        assert np.all(
            ion_solution >= 0.0
        ), "Negative ion number density found, this should not happen."
        assert (
            electron_solution >= 0.0
        ), "Negative electron density found, this should not happen."
        delta_ion = (ion_number_density - ion_solution) / ion_solution
        delta_electron = (
            electron_densities - electron_solution
        ) / electron_solution
        return delta_ion, delta_electron

    @staticmethod
    def solve(rate_matrix, balance_vector):
        """Solves the NLTE ionization rate equations.

        Parameters
        ----------
        rate_matrix : pandas.DataFrame
            Rate matrix used for the solver.
        balance_vector : numpy.array
            Balance vector used for the solver.

        Returns
        -------
        numpy.array
            Solution vector for the NLTE ionization solver.
        """
        lu, piv = lu_factor(rate_matrix)
        return lu_solve((lu, piv), balance_vector)

    @staticmethod
    def calculate_balance_vector(
        number_density,
        rate_matrix_index,
    ):
        """Constructs the balance vector for the NLTE ionization solver set of equations by combining
        all solution verctor blocks.

        Parameters
        ----------
        number_density : pandas.DataFrame
            Number densities of all present species.
        rate_matrix_index : pandas.MultiIndex
            (atomic_number, ion_number, treatment type)

        Returns
        -------
        numpy.array
            Solution vector for the NLTE ionization solver.
        """
        atomic_numbers = number_density.index
        balance_array = []
        for atomic_number in atomic_numbers:
            needed_number_of_elements = (
                rate_matrix_index.get_level_values("atomic_number")
                == atomic_number
            ).sum()
            balance_vector_block = np.zeros(needed_number_of_elements)
            balance_vector_block[-1] = number_density.loc[atomic_number]
            balance_array.append(balance_vector_block)
        return np.hstack(balance_array)

    @staticmethod
    def prepare_ion_recomb_coefficients_nlte_ion(
        gamma,
        alpha_sp,
        alpha_stim,
        coll_ion_coeff,
        coll_recomb_coeff,
        partition_function,
        levels,
        level_boltzmann_factor,
    ):
        """
        Prepares the ionization and recombination coefficients by grouping them for
        ion numbers.

        Parameters
        ----------
        gamma : pandas.DataFrame
            The rate coefficient for radiative ionization.
        alpha_sp : pandas.DataFrame
            The rate coefficient for spontaneous recombination.
        alpha_stim : pandas.DataFrame
            The rate coefficient for stimulated recombination.
        coll_ion_coeff : pandas.DataFrame
            The rate coefficient for collisional ionization in the Seaton
            approximation.
        coll_recomb_coeff : pandas.DataFrame
            The rate coefficient for collisional recombination.
        partition_function : pandas.DataFrame
            General partition function. Indexed by atomic number, ion number.
        levels : MultiIndex
            (atomic_number, ion_number, level_number)
            Index of filtered atomic data.
        level_boltzmann_factor : pandas.DataFrame
            General Boltzmann factor.
        Returns
        -------
        total_photo_ion_coefficients
            Photoionization coefficients grouped by atomic number and ion number.
        total_rad_recomb_coefficients
            Radiative recombination coefficients grouped by atomic number and ion number.
        total_coll_ion_coefficients
            Collisional ionization coefficients grouped by atomic number and ion number.
        total_coll_recomb_coefficients
            Collisional recombination coefficients grouped by atomic number and ion number.
        """
        indexer = pd.Series(
            np.arange(partition_function.shape[0]),
            index=partition_function.index,
        )
        _ion2level_idx = indexer.loc[levels.droplevel("level_number")].values
        partition_function_broadcast = partition_function.values[_ion2level_idx]
        level_population_fraction = pd.DataFrame(
            level_boltzmann_factor.values / partition_function_broadcast,
            index=levels,
        )
        total_photo_ion_coefficients = (
            (level_population_fraction.loc[gamma.index] * gamma)
            .groupby(level=("atomic_number", "ion_number"))
            .sum()
        )

        total_rad_recomb_coefficients = (
            (alpha_sp + alpha_stim)
            .groupby(level=["atomic_number", "ion_number"])
            .sum()
        )
        total_coll_ion_coefficients = (
            (
                level_population_fraction.loc[coll_ion_coeff.index]
                * coll_ion_coeff
            )
            .groupby(level=("atomic_number", "ion_number"))
            .sum()
        )
        total_coll_recomb_coefficients = (
            (coll_recomb_coeff)
            .groupby(level=("atomic_number", "ion_number"))
            .sum()
        )
        return (
            total_photo_ion_coefficients,
            total_rad_recomb_coefficients,
            total_coll_ion_coefficients,
            total_coll_recomb_coefficients,
        )

    @staticmethod
    def calculate_rate_matrix(
        atomic_numbers,
        phi_shell,
        electron_density,
        rate_matrix_index,
        total_photo_ion_coefficients,
        total_rad_recomb_coefficients,
        total_coll_ion_coefficients,
        total_coll_recomb_coefficients,
    ):
        """

        Parameters
        ----------
        phi_shell : pandas.DataFrame
            Saha Factors in the current shell
        electron_density : float
            Guess for electron density in the current shell
        rate_matrix_index : pandas.MultiIndex
            Index used for constructing the rate matrix
        total_photo_ion_coefficients : pandas.DataFrame
            Photo ionization coefficients
        total_rad_recomb_coefficients : pandas.DataFrame
            Radiative recombination coefficients (should get multiplied by electron density)
        total_coll_ion_coefficients : pandas.DataFrame
            Collisional ionization coefficients (should get multiplied by electron density)
        total_coll_recomb_coefficients : pandas.DataFrame
            Collisional recombination coefficients (should get multiplied by electron density^2)

        Returns
        -------
        pandas.DataFrame
            Rate matrix used for NLTE solver.
        """
        rate_matrix = pd.DataFrame(
            0.0, columns=rate_matrix_index, index=rate_matrix_index
        )
        total_rad_recomb_coefficients = (
            total_rad_recomb_coefficients * electron_density
        )
        total_coll_ion_coefficients = (
            total_coll_ion_coefficients * electron_density
        )
        total_coll_recomb_coefficients = (
            total_coll_recomb_coefficients * electron_density**2
        )
        for atomic_number in atomic_numbers:
            ion_numbers = rate_matrix.loc[atomic_number].index.get_level_values(
                "ion_number"
            )
            phi_block = phi_shell.loc[atomic_number]
            rate_matrix_block = NLTEPopulationSolver.lte_rate_matrix_block(
                phi_block, electron_density
            )

            nlte_ion_numbers = ion_numbers[
                rate_matrix.loc[atomic_number].index.get_level_values(
                    "level_number"
                )
                == "nlte_ion"
            ]
            # >>> lte_ion_numbers is for future use in NLTE excitation treatment
            lte_ion_numbers = ion_numbers[
                rate_matrix.loc[atomic_number].index.get_level_values(
                    "level_number"
                )
                == "lte_ion"
            ]
            # <<<
            for ion_number in nlte_ion_numbers:
                rate_matrix_block = NLTEPopulationSolver.set_nlte_ion_rate(
                    rate_matrix_block,
                    atomic_number,
                    ion_number,
                    total_rad_recomb_coefficients.loc[(atomic_number,)],
                    total_photo_ion_coefficients.loc[(atomic_number,)],
                    total_coll_ion_coefficients.loc[(atomic_number,)],
                    total_coll_recomb_coefficients.loc[(atomic_number,)],
                )
            rate_matrix.loc[
                (atomic_number, slice(None)), (atomic_number)
            ] = rate_matrix_block
        return rate_matrix

    @staticmethod
    def set_nlte_ion_rate(
        rate_matrix_block,
        atomic_number,
        ion_number,
        total_rad_recomb_coefficients,
        total_photo_ion_coefficients,
        total_coll_ion_coefficients,
        total_coll_recomb_coefficients,
    ):
        """Calculates the row for the species treated in NLTE ionization

        Parameters
        ----------
        rate_matrix_block : numpy.array
            The diagonal block corresponding to current atomic number.
        atomic_number : int
            Current atomic number
        ion_number : int
            Current ion number
        total_rad_recomb_coefficients : pandas.DataFrame
            Rad. recomb. coefficients for current atomic number
        total_photo_ion_coefficients : pandas.DataFrame
            Photo ionization coefficients for current atomic number
        total_coll_ion_coefficients : pandas.DataFrame
            Collisional ionization coefficients for current atomic number
        total_coll_recomb_coefficients : pandas.DataFrame
            Collisional recombination coefficients for current atomic number

        Returns
        -------
        numpy.array
            Rate matrix block with a changed row for NLTE ionization treatment
        """
        ion_coefficients = (
            total_photo_ion_coefficients + total_coll_ion_coefficients
        )
        recomb_coefficients = (
            total_rad_recomb_coefficients + total_coll_recomb_coefficients
        )
        if atomic_number != ion_number:
            ion_coeff_matrix_ion_row = NLTEPopulationSolver.ion_matrix(
                ion_coefficients, atomic_number, ion_number
            )
            recomb_coeff_matrix_ion_row = NLTEPopulationSolver.recomb_matrix(
                recomb_coefficients, atomic_number, ion_number
            )
            rate_matrix_block[ion_number, :] = (
                ion_coeff_matrix_ion_row + recomb_coeff_matrix_ion_row
            )
        return rate_matrix_block

    @staticmethod
    def lte_rate_matrix_block(phi_block, electron_density):
        """Creates the generic LTE block for rate matrix.

        Parameters
        ----------
        phi_block : pandas.DataFrame
            Saha Factors for current atomic number
        electron_density : float
            Current guess for electron density

        Returns
        -------
        numpy.array
            LTE block for rate matrix
        """
        lte_rate_vector_block = -1.0 * np.hstack([*phi_block.values, -1.0])
        lte_rate_matrix_block = np.diag(lte_rate_vector_block)
        n_e_initial = np.ones(len(phi_block)) * electron_density
        n_e_matrix = np.diag(n_e_initial, 1)
        lte_rate_matrix_block += n_e_matrix
        lte_rate_matrix_block[-1, :] = 1.0
        return lte_rate_matrix_block

    @staticmethod
    def recomb_matrix(recomb_coefficients, atomic_number, ion_number):
        """Constructs a recombination rate matrix from the recombination rates.

        Parameters
        ----------
        recomb_coefficients : pandas.DataFrame
            Recombination coefficients.
        atomic_number : int64
            Current atomic number. Used for the dimension of a square matrix.
        ion_number : int64
            Current ion number. Used for returning the correct row.

        Returns
        -------
        numpy.ndarray
        """
        offdiag = np.zeros(atomic_number)
        index = recomb_coefficients.index
        for i in index:
            offdiag[i] = recomb_coefficients.loc[i]
        diag = np.hstack([np.zeros(1), -offdiag])
        return (np.diag(diag) + np.diag(offdiag, k=1))[ion_number, :]

    @staticmethod
    def ion_matrix(ion_coefficients, atomic_number, ion_number):
        """Constructs an ionization rate matrix from the ionization rates.

        Parameters
        ----------
        ion_coefficients : pandas.DataFrame
            Recombination coefficients.
        atomic_number : int64
            Current atomic number. Used for the dimension of a square matrix.
        ion_number : int64
            Current ion number. Used for returning the correct row.

        Returns
        -------
        numpy.ndarray
        """
        offdiag = np.zeros(atomic_number)
        index = ion_coefficients.index
        for i in index:
            offdiag[i] = ion_coefficients.loc[i]
        diag = np.hstack([-offdiag, np.zeros(1)])
        return (np.diag(diag) + np.diag(offdiag, k=-1))[ion_number, :]
