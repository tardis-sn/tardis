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
        alpha_sp : DataFrame
        alpha_stim : DataFrame
        coll_ion_coeff : DataFrame
        coll_recomb_coeff : DataFrame
        partition_function : DataFrame
        levels : MultiIndex
        level_boltzmann_factor : DataFrame
        phi : DataFrame
        rate_matrix_index : MultiIndex
        number_density : DataFrame

        Returns
        -------
        ion_number_densities_nlte: DataFrame
        electron_densities_nlte: Series
        """

        (
            photo_ion_rate,
            rad_recomb_rate_coeff,
            coll_ion_coefficient,
            coll_recomb_coefficient,
        ) = NLTERateEquationSolver.prepare_ion_recomb_rates_nlte_ion(
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
        rate_matrix = self.calculate_rate_matrix(
            phi[0],
            initial_electron_density[0],
            rate_matrix_index,
            photo_ion_rate[0],
            rad_recomb_rate_coeff[0],
            coll_ion_coefficient[0],
            coll_recomb_coefficient[0],
        )
        return (
            -1,
            -1,
        )  # function still empty, that's why return statement is arbitrary at this point

    @staticmethod
    def calculate_rate_matrix(
        phi_shell,
        electron_density,
        rate_matrix_index,
        photo_ion_rate,
        rad_recomb_rate_coeff,
        coll_ion_coefficient,
        coll_recomb_coefficient,
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
        photo_ion_rate : DataFrame
            Photo ionization rates
        rad_recomb_rate_coeff : DataFrame
            Radiative recombination coefficients(should get multiplied by electron density)
        coll_ion_coefficient : DataFrame
            Collisionional ionization coefficients(should get multiplied by electron density)
        coll_recomb_coefficient : DataFrame
            Collisional recombination coefficients (should get multiplied by electron density^2)

        Returns
        -------
        DataFrame
            Rate matrix used for nlte solver.
        """
        rate_matrix = pd.DataFrame(
            0.0, columns=rate_matrix_index, index=rate_matrix_index
        )
        rad_recomb_rates = rad_recomb_rate_coeff * electron_density
        coll_ion_rates = coll_ion_coefficient * electron_density
        coll_recomb_rates = coll_recomb_coefficient * electron_density**2
        atomic_numbers = (
            rate_matrix_index.get_level_values(0).unique().drop("n_e")
        )  # dropping the n_e index
        for atomic_number in atomic_numbers:
            ion_numbers = rate_matrix.loc[atomic_number].index.get_level_values(
                0
            )
            phi_block = phi_shell.loc[atomic_number]
            rate_matrix_block = NLTERateEquationSolver.lte_rate_matrix_block(
                phi_block, electron_density
            )

            nlte_ion_numbers = ion_numbers[
                rate_matrix.loc[atomic_number].index.get_level_values(1)
                == "nlte_ion"
            ]
            lte_ion_numbers = ion_numbers[
                rate_matrix.loc[atomic_number].index.get_level_values(1)
                == "lte_ion"
            ]
            for ion_number in nlte_ion_numbers:
                rate_matrix_block = NLTERateEquationSolver.set_nlte_ion_rate(
                    rate_matrix_block,
                    atomic_number,
                    ion_number,
                    rad_recomb_rates.loc[(atomic_number,)],
                    photo_ion_rate.loc[(atomic_number,)],
                    coll_ion_rates.loc[(atomic_number,)],
                    coll_recomb_rates.loc[(atomic_number,)],
                )
            rate_matrix.loc[
                (atomic_number, slice(None)), (atomic_number)
            ] = rate_matrix_block

        last_row = NLTERateEquationSolver.prepare_last_row(atomic_numbers)
        rate_matrix.loc[("n_e", slice(None))] = last_row
        return rate_matrix

    @staticmethod
    def set_nlte_ion_rate(
        rate_matrix_block,
        atomic_number,
        ion_number,
        radiative_recombination_rate,
        photo_ion_rates,
        coll_ion_rate,
        coll_recomb_rate,
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
        radiative_recombination_rate : DataFrame
            Rad. recomb. rate for current atomic number
        photo_ion_rates : DataFrame
            Photo ion. rate for current atomic number
        coll_ion_rate : DataFrame
            Coll. ion. rate for current atomic number
        coll_recomb_rate : DataFrame
            Coll. recomb. rate for current atomic number

        Returns
        -------
        numpy.array
            Rate matrix block with a changed row for NLTE ionization treatment
        """
        ion_rates = photo_ion_rates + coll_ion_rate
        recomb_rate = radiative_recombination_rate + coll_recomb_rate
        if atomic_number != ion_number:
            ion_rate_matrix = NLTERateEquationSolver.ion_matrix(
                ion_rates, atomic_number
            )
            recomb_rate_matrix = NLTERateEquationSolver.recomb_matrix(
                recomb_rate, atomic_number
            )
            rate_matrix_block[ion_number, :] = (
                ion_rate_matrix + recomb_rate_matrix
            )[ion_number, :]
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
    def recomb_matrix(recomb_rate, atomic_number):
        """Constructs a recombination rate matrix from the recombination rates.

        Parameters
        ----------
        recomb_rate : DataFrame
            Recombination rates.
        atomic_number : int64
            Current atomic number. Used for the dimension of a square matrix.

        Returns
        -------
        numpy.ndarray
        """
        offdiag = np.zeros(atomic_number)
        index = recomb_rate.index
        for i in index:
            offdiag[i] = recomb_rate.loc[i]
        diag = np.hstack([np.zeros(1), -offdiag])
        return np.diag(diag) + np.diag(offdiag, k=1)

    @staticmethod
    def ion_matrix(ion_rate, atomic_number):
        """Constructs an ionization rate matrix from the ionization rates.

        Parameters
        ----------
        recomb_rate : DataFrame
            Recombination rates.
        atomic_number : int64
            Current atomic number. Used for the dimension of a square matrix.

        Returns
        -------
        numpy.ndarray
        """
        offdiag = np.zeros(atomic_number)
        index = ion_rate.index
        for i in index:
            offdiag[i] = ion_rate.loc[i]
        diag = np.hstack([-offdiag, np.zeros(1)])
        return np.diag(diag) + np.diag(offdiag, k=-1)

    @staticmethod
    def prepare_last_row(atomic_numbers):
        """Prepares the last row of the rate_matrix. This row corresponds to the charge density equation."""
        last_row = []
        for atomic_number in atomic_numbers:
            last_row.append(np.arange(0.0, atomic_number + 1))
        last_row = np.hstack([*last_row, -1])
        # TODO needs to be modified for use in nlte_excitation
        return last_row

    @staticmethod
    def prepare_ion_recomb_rates_nlte_ion(
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
        Prepares the ionization and recombination rates/coefficients by grouping them for ion numbers.
        """
        indexer = pd.Series(
            np.arange(partition_function.shape[0]),
            index=partition_function.index,
        )
        _ion2level_idx = indexer.loc[levels.droplevel(2)].values
        partition_function_broadcast = partition_function.values[_ion2level_idx]
        level_population_fraction = pd.DataFrame(
            level_boltzmann_factor.values / partition_function_broadcast,
            index=levels,
        )
        photo_ion_rate = (
            (level_population_fraction.loc[gamma.index] * gamma)
            .groupby(level=(0, 1))
            .sum()
        )
        rad_recomb_rate_coeff = (
            alpha_sp.groupby(level=[0, 1]).sum()
            + alpha_stim.groupby(level=[0, 1]).sum()
        )
        coll_ion_coefficient = (
            (
                level_population_fraction.loc[coll_ion_coeff.index]
                * coll_ion_coeff
            )
            .groupby(level=(0, 1))
            .sum()
        )
        coll_recomb_coefficient = (
            (coll_recomb_coeff).groupby(level=(0, 1)).sum()
        )
        return (
            photo_ion_rate,
            rad_recomb_rate_coeff,
            coll_ion_coefficient,
            coll_recomb_coefficient,
        )
