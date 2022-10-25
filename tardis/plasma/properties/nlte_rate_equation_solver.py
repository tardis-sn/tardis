import pandas as pd
import numpy as np

from tardis.plasma.properties.base import ProcessingPlasmaProperty

__all__ = [
    "NLTERateEquationSolver",
]


class NLTERateEquationSolver(ProcessingPlasmaProperty):
    outputs = ("ion_number_density_nlte", "n_e")
    def calculate(self, *args, **kwargs):
        return super().calculate(*args, **kwargs)


    def calculate_rate_matrix(
        self,
        phi,
        electron_density,
        number_density,
        rate_matrix_index,
        last_row,
        rad_recomb_rate_coeff,
        photo_ion_rates,
        coll_ion_coefficient,
        coll_recomb_coefficient,
    ):
        1/0
        rate_matrix = pd.DataFrame(
            0.0, columns=rate_matrix_index, index=rate_matrix_index
        )
        rad_recomb_rates = rad_recomb_rate_coeff * electron_density
        coll_ion_rates = coll_ion_coefficient * electron_density
        coll_recomb_rates= coll_recomb_coefficient * electron_density**2
        atomic_numbers = number_density.index
        for atomic_number in atomic_numbers:
            ion_numbers = rate_matrix.loc[atomic_number].index.get_level_values(0)
            phi_block = phi.loc[atomic_number]
            rate_matrix_block = self.lte_rate_matrix_block(phi_block, electron_density)

            nlte_ion_numbers = ion_numbers[
                rate_matrix.loc[atomic_number].index.get_level_values(1) == "nlte_ion"
            ]
            lte_ion_numbers = ion_numbers[
                rate_matrix.loc[atomic_number].index.get_level_values(1) == "lte_ion"
            ]
            for ion_number in nlte_ion_numbers:
                rate_matrix_block = self.set_nlte_ion_rate(
                    rate_matrix_block,
                    atomic_number,
                    ion_number,
                    rad_recomb_rates.loc[(atomic_number,)],
                    photo_ion_rates.loc[(atomic_number,)],
                    coll_ion_rates.loc[(atomic_number,)],
                    coll_recomb_rates.loc[(atomic_number,)],
                )
            rate_matrix.loc[
                (atomic_number, slice(None)), (atomic_number)
            ] = rate_matrix_block

        rate_matrix.loc[("n_e", slice(None))] = last_row
        return rate_matrix


    def set_nlte_ion_rate(
        self,
        rate_matrix_block,
        atomic_number,
        ion_number,
        radiative_recombination_rate,
        photo_ion_rates,
        coll_ion_rate,
        coll_recomb_rate,
    ):
        ion_rates = photo_ion_rates + coll_ion_rate
        recomb_rate = radiative_recombination_rate + coll_recomb_rate
        if atomic_number == ion_number:
            rate_matrix_block[ion_number, :] = 1.0
        else:
            ion_rate_matrix = self.ion_matrix(ion_rates, atomic_number)
            recomb_rate_matrix = self.recomb_matrix(recomb_rate, atomic_number)
            rate_matrix_block[ion_number, :] = (ion_rate_matrix + recomb_rate_matrix)[
                ion_number, :
            ]
        return rate_matrix_block


    def lte_rate_matrix_block(self, phi_block, electron_density):
        lte_rate_vector_block = np.hstack([-phi_block, 1.0])
        lte_rate_matrix_block = np.diag(lte_rate_vector_block)
        n_e_initial = np.ones(len(phi_block)) * electron_density
        n_e_matrix = np.diag(n_e_initial, 1)
        lte_rate_matrix_block += n_e_matrix
        lte_rate_matrix_block[-1, :] = 1.0
        return lte_rate_matrix_block


    def prepare_phi(self, phi):
        phi[phi == 0.0] = 1.0e-10 * phi[phi > 0.0].min().min()
        return phi


    def recomb_matrix(self, recomb_rate, atomic_number):
        offdiag = np.zeros(atomic_number)
        index = recomb_rate.index
        for i in index:
            offdiag[i] = recomb_rate[i]
        diag = np.hstack([np.zeros(1), -offdiag])
        return np.diag(diag) + np.diag(offdiag, k=1)


    def ion_matrix(self, ion_rate, atomic_number):
        offdiag = np.zeros(atomic_number)
        index = ion_rate.index
        for i in index:
            offdiag[i] = ion_rate[i]
        diag = np.hstack([-offdiag, np.zeros(1)])
        return np.diag(diag) + np.diag(offdiag, k=-1)


    @staticmethod
    def prepare_last_row(atomic_numbers):
        last_row = []
        for atomic_number in atomic_numbers:
            last_row.append(np.arange(0.0, atomic_number + 1))
        last_row = np.hstack([*last_row, -1])
        return last_row


    def prepare_ion_recomb_rates_nlte_ion(
        self,
        level_population_fraction,
        gamma,
        alpha_sp,
        alpha_stim,
        coll_ion_coeff,
        coll_recomb_coeff,
        partition_function,
        levels,
        level_boltzmann_factor,
    ):
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
            (level_population_fraction.loc[gamma.index] * gamma).groupby(level=(0, 1)).sum()
        )
        rad_recomb_rate_coeff = (
            alpha_sp.groupby(level=[0, 1]).sum() + alpha_stim.groupby(level=[0, 1]).sum()
        )
        coll_ion_coefficient = (
            (level_population_fraction.loc[coll_ion_coeff.index] * coll_ion_coeff)
            .groupby(level=(0, 1))
            .sum()
        )
        coll_recomb_coefficient = (coll_recomb_coeff).groupby(level=(0, 1)).sum()
        return (
            photo_ion_rate,
            rad_recomb_rate_coeff,
            coll_ion_coefficient,
            coll_recomb_coefficient,
        )
