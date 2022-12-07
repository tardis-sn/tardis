import pandas as pd
import numpy as np

from tardis.plasma.properties.base import ProcessingPlasmaProperty

__all__ = [
    "NLTERateEquationSolver",
]


class NLTERateEquationSolver(ProcessingPlasmaProperty):
    outputs = ("ion_number_density_nlte", "electron_densities_nlte")

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
    ):
        """Calculates ion number densities and electron densities using NLTE ionization.

        Parameters
        ----------
        gamma : DataFrame
            The rate coefficient for radiative ionization.
        alpha_sp : DataFrame
            The rate coefficient for spontaneous recombination.
        alpha_stim : DataFrame
            The rate coefficient for stimulated recombination.
        coll_ion_coeff : DataFrame
            The rate coefficient for collisional ionization in the Seaton
            approximation.
        coll_recomb_coeff : DataFrame
            The rate coefficient for collisional recombination.
        partition_function : DataFrame
            General partition function. Indexed by atomic number, ion number.
        levels : MultiIndex
            (atomic_number, ion_number, level_number)
            Index of filtered atomic data.
        level_boltzmann_factor : DataFrame
            General Boltzmann factor.
        phi : DataFrame
            Saha Factors.
        rate_matrix_index : MultiIndex
            (atomic_number, ion_number, treatment type)
            If ion is treated in LTE or nebular ionization, 3rd index is "lte_ion",
            if treated in NLTE ionization, 3rd index is "nlte_ion".
        number_density : DataFrame
            Number density in each shell for each species.

        Returns
        -------
        ion_number_densities_nlte : DataFrame
            Number density with NLTE ionization treatment.
        electron_densities_nlte : Series
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

        # >>>TODO:initial electron density should be included in the initial guess, added in a future PR
        initial_electron_density = number_density.sum(axis=0)
        # <<<
        atomic_numbers = (
            rate_matrix_index.get_level_values("atomic_number")
            .unique()
            .drop("n_e")
        )  # dropping the n_e index, as rate_matrix_index's first index is (atomic_numbers, "n_e")
        rate_matrix = self.calculate_rate_matrix(
            atomic_numbers,
            phi[0],
            initial_electron_density[0],
            rate_matrix_index,
            total_photo_ion_coefficients[0],
            total_rad_recomb_coefficients[0],
            total_coll_ion_coefficients[0],
            total_coll_recomb_coefficients[0],
        )
        initial_guess = self.prepare_first_guess(
            atomic_numbers, number_density[0], initial_electron_density[0]
        )
        jacobian_matrix = self.jacobian_matrix(
            atomic_numbers,
            initial_guess,
            rate_matrix,
            rate_matrix_index,
            total_rad_recomb_coefficients[0],
            total_coll_ion_coefficients[0],
            total_coll_recomb_coefficients[0],
        )
        # TODO: change the jacobian and rate matrix to use shell id and get coefficients from the attribute of the class.

        raise NotImplementedError(
            "NLTE ionization hasn't been fully implemented yet!"
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
        phi_shell : DataFrame
            Saha Factors in the current shell
        electron_density : float
            Guess for electron density in the current shell
        rate_matrix_index : MultiIndex
            Index used for constructing the rate matrix
        total_photo_ion_coefficients : DataFrame
            Photo ionization coefficients
        total_rad_recomb_coefficients : DataFrame
            Radiative recombination coefficients (should get multiplied by electron density)
        total_coll_ion_coefficients : DataFrame
            Collisional ionization coefficients (should get multiplied by electron density)
        total_coll_recomb_coefficients : DataFrame
            Collisional recombination coefficients (should get multiplied by electron density^2)

        Returns
        -------
        DataFrame
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
            rate_matrix_block = NLTERateEquationSolver.lte_rate_matrix_block(
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
                rate_matrix_block = NLTERateEquationSolver.set_nlte_ion_rate(
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

        charge_conservation_row = (
            NLTERateEquationSolver.prepare_charge_conservation_row(
                atomic_numbers
            )
        )
        rate_matrix.loc[("n_e", slice(None))] = charge_conservation_row
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
        total_rad_recomb_coefficients : DataFrame
            Rad. recomb. coefficients for current atomic number
        total_photo_ion_coefficients : DataFrame
            Photo ion. coefficients for current atomic number
        total_coll_ion_coefficients : DataFrame
            Coll. ion. coefficients for current atomic number
        total_coll_recomb_coefficients : DataFrame
            Coll. recomb. coefficients for current atomic number

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
            ion_coeff_matrix_ion_row = NLTERateEquationSolver.ion_matrix(
                ion_coefficients, atomic_number, ion_number
            )
            recomb_coeff_matrix_ion_row = NLTERateEquationSolver.recomb_matrix(
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
        phi_block : DataFrame
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
    def prepare_phi(phi):
        """
        Makes sure that phi does not have any 0 entries.
        """
        phi[phi == 0.0] = 1.0e-10 * phi[phi > 0.0].min().min()
        return phi

    @staticmethod
    def recomb_matrix(recomb_coefficients, atomic_number, ion_number):
        """Constructs a recombination rate matrix from the recombination rates.

        Parameters
        ----------
        recomb_coefficients : DataFrame
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
        ion_coefficients : DataFrame
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

    @staticmethod
    def prepare_charge_conservation_row(atomic_numbers):
        """Prepares the last row of the rate_matrix. This row corresponds to the charge
        density equation."""
        charge_conservation_row = []
        for atomic_number in atomic_numbers:
            charge_conservation_row.append(np.arange(0, atomic_number + 1))
        charge_conservation_row = np.hstack([*charge_conservation_row, -1])
        # TODO needs to be modified for use in nlte_excitation
        return charge_conservation_row

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
        gamma : DataFrame
            The rate coefficient for radiative ionization.
        alpha_sp : DataFrame
            The rate coefficient for spontaneous recombination.
        alpha_stim : DataFrame
            The rate coefficient for stimulated recombination.
        coll_ion_coeff : DataFrame
            The rate coefficient for collisional ionization in the Seaton
            approximation.
        coll_recomb_coeff : DataFrame
            The rate coefficient for collisional recombination.
        partition_function : DataFrame
            General partition function. Indexed by atomic number, ion number.
        levels : MultiIndex
            (atomic_number, ion_number, level_number)
            Index of filtered atomic data.
        level_boltzmann_factor : DataFrame
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
    def jacobian_matrix(
        atomic_numbers,
        populations,
        rate_matrix,
        rate_matrix_index,
        total_rad_recomb_coefficients,
        total_coll_ion_coefficients,
        total_coll_recomb_coefficients,
    ):
        """Creates the jacobian matrix used for NLTE ionization solver

        Parameters
        ----------
        populations : numpy.array
            Ion populations, electron density
        rate_matrix : DataFrame
            Rate matrix used for NLTE solver.
        rate_matrix_index : MultiIndex
            (atomic_number, ion_number, treatment type)
            If ion is treated in LTE or nebular ionization, 3rd index is "lte_ion",
            if treated in NLTE ionization, 3rd index is "nlte_ion".
        total_rad_recomb_coefficients : DataFrame
            Radiative recombination coefficients grouped by atomic number and ion number.
        total_coll_ion_coefficients : DataFrame
            Collisional ionization coefficients(should get multiplied by electron density).
        total_coll_recomb_coefficients : DataFrame
            Collisional recombination coefficients(should get multiplied by electron density).

        Returns
        -------
        numpy.array
            Jacobian matrix used for NLTE ionization solver
        """
        # TODO: for future use, can be vectorized.
        index = 0
        jacobian_matrix = rate_matrix.copy().values
        jacobian_matrix[:-1, -1] = populations[1:]
        for atomic_number in atomic_numbers:
            for i in range(index, index + atomic_number):
                if rate_matrix_index[i][2] == "nlte_ion":
                    jacobian_matrix[
                        i, -1
                    ] = NLTERateEquationSolver.deriv_matrix_block(
                        atomic_number,
                        total_rad_recomb_coefficients.loc[(atomic_number,)],
                        total_coll_ion_coefficients.loc[(atomic_number,)],
                        total_coll_recomb_coefficients.loc[(atomic_number,)],
                        populations[index : index + atomic_number + 1],
                        populations[-1],
                    )[
                        i - index
                    ]
            index += atomic_number + 1
            jacobian_matrix[index - 1, -1] = 0  # number conservation row
        return jacobian_matrix

    @staticmethod
    def deriv_matrix_block(
        atomic_number,
        total_rad_recomb_coefficients,
        total_coll_ion_coefficients,
        total_coll_recomb_coefficients,
        current_ion_number_densities,
        current_electron_density,
    ):
        """Calculates the dot product of the derivative of rate matrix and ion number densities+electron density column.

        Parameters
        ----------
        atomic_number : int64
            Current atomic number
        total_rad_recomb_coefficients : DataFrame
            Radiative recombination coefficients grouped by atomic number and ion number.
        total_coll_ion_coefficients : DataFrame
            Collisional ionization coefficients.
        total_coll_recomb_coefficients : DataFrame
            Collisional recombination coefficients.
        current_ion_number_densities : numpy.array
            Current ion number densities for the current atomic number.
        current_electron_density : float64
            Current electron density

        Returns
        -------
        numpy.array
            Returns the part of the last column of the jacobian matrix, corresponding to atomic number.
        """
        ion_numbers = np.arange(0, atomic_number)
        radiative_rate_coeff_matrix = NLTERateEquationSolver.recomb_matrix(
            total_rad_recomb_coefficients, atomic_number, ion_numbers
        )
        coll_recomb_matrix = (
            NLTERateEquationSolver.recomb_matrix(
                total_coll_recomb_coefficients, atomic_number, ion_numbers
            )
            * current_electron_density
            * 2
        )
        coll_ion_coeff_matrix = NLTERateEquationSolver.ion_matrix(
            total_coll_ion_coefficients, atomic_number, ion_numbers
        )
        deriv_matrix = (
            radiative_rate_coeff_matrix
            + coll_ion_coeff_matrix
            + coll_recomb_matrix
        )
        return np.dot(deriv_matrix, current_ion_number_densities)

    def prepare_first_guess(
        self, atomic_numbers, number_density, electron_density
    ):
        # TODO needs to be changed for excitation
        array_size = (number_density.index.values + 1).sum() + 1
        first_guess = np.zeros(array_size)
        index = 1
        for atomic_number in atomic_numbers:
            first_guess[index] = number_density.loc[atomic_number]
            index += atomic_number + 1
        first_guess[-1] = electron_density
        return first_guess
