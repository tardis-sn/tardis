import logging
import warnings

import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.optimize import root

from tardis.iip_plasma.continuum.util import *
from tardis.iip_plasma.exceptions import PlasmaIonizationError
from tardis.iip_plasma.properties.base import ProcessingPlasmaProperty
from tardis.iip_plasma.properties.partition_function import PartitionFunction

logger = logging.getLogger(__name__)

__all__ = [
    "IonNumberDensity",
    "LTEIonNumberDensity",
    "NLTEIonNumberDensity",
    "PhiSahaElectrons",
    "PhiSahaLTE",
    "PhiSahaLTECont",
    "PhiSahaNebular",
    "RadiationFieldCorrection",
]


def calculate_block_ids_from_dataframe(dataframe):
    block_start_id = (
        np.where(np.diff(dataframe.index.get_level_values(0)) != 0.0)[0] + 1
    )
    return np.hstack(([0], block_start_id, [len(dataframe)]))


class PhiSahaLTE(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    phi : Pandas DataFrame, dtype float
          Used for LTE ionization. Indexed by atomic number, ion number. Columns are zones.
    """

    outputs = ("phi",)
    latex_name = ("\\Phi",)
    latex_formula = (
        "\\dfrac{2Z_{i,j+1}}{Z_{i,j}}\\Big(\
                     \\dfrac{2\\pi m_{e}/\\beta_{\\textrm{rad}}}{h^2}\
                     \\Big)^{3/2}e^{\\dfrac{-\\chi_{i,j}}{kT_{\
                     \\textrm{rad}}}}",
    )

    broadcast_ionization_energy = None

    @staticmethod
    def calculate(g_electron, beta_rad, partition_function, ionization_data):
        # print "Calculating Phi LTE"

        phis = np.empty(
            (
                partition_function.shape[0]
                - partition_function.index.get_level_values(0).unique().size,
                partition_function.shape[1],
            )
        )

        block_ids = calculate_block_ids_from_dataframe(partition_function)

        for i, start_id in enumerate(block_ids[:-1]):
            end_id = block_ids[i + 1]
            current_block = partition_function.values[start_id:end_id]
            current_phis = current_block[1:] / current_block[:-1]
            phis[start_id - i : end_id - i - 1] = current_phis

        broadcast_ionization_energy = ionization_data.reindex(
            partition_function.index
        ).dropna()
        phi_index = broadcast_ionization_energy.index
        broadcast_ionization_energy = broadcast_ionization_energy.values

        phi_coefficient = (
            2
            * g_electron
            * np.exp(np.outer(broadcast_ionization_energy, -beta_rad))
        )

        return pd.DataFrame(phis * phi_coefficient, index=phi_index).fillna(0.0)

    @staticmethod
    def _calculate_block_ids(partition_function):
        partition_function.index.get_level_values(0).unique()


class PhiSahaLTECont(PhiSahaLTE):
    outputs = ("phi_lte",)
    latex_name = ("\\Phi_lte",)


class PhiSahaElectrons(PhiSahaLTE):
    outputs = ("phi_Te",)
    latex_name = ("\\Phi_Te",)

    @staticmethod
    def calculate(
        g_electron_Te, beta_electron, lte_partition_function_Te, ionization_data
    ):
        return super(PhiSahaElectrons, PhiSahaElectrons).calculate(
            g_electron_Te,
            beta_electron,
            lte_partition_function_Te,
            ionization_data,
        )


class PhiSahaNebular(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    phi : Pandas DataFrame, dtype float
          Used for nebular ionization. Indexed by atomic number, ion number. Columns are zones.
    """

    outputs = ("phi",)
    latex_name = ("\\Phi",)
    latex_formula = (
        "W(\\delta\\zeta_{i,j}+W(1-\\zeta_{i,j}))\\left(\
                     \\dfrac{T_{\\textrm{electron}}}{T_{\\textrm{rad}}}\
                     \\right)^{1/2}",
    )

    @staticmethod
    def calculate(
        t_rad,
        w,
        zeta_data,
        t_electrons,
        delta,
        g_electron,
        beta_rad,
        partition_function,
        ionization_data,
    ):
        # print "Calculating Phi Nebular"
        phi_lte = PhiSahaLTE.calculate(
            g_electron, beta_rad, partition_function, ionization_data
        )
        zeta = PhiSahaNebular.get_zeta_values(zeta_data, phi_lte.index, t_rad)
        phis = (
            phi_lte
            * w
            * ((zeta * delta) + w * (1 - zeta))
            * (t_electrons / t_rad) ** 0.5
        )
        return phis.fillna(0.0)

    @staticmethod
    def get_zeta_values(zeta_data, ion_index, t_rad):
        zeta_t_rad = zeta_data.columns.values.astype(np.float64)
        zeta_values = zeta_data.loc[ion_index].values.astype(np.float64)
        zeta = interpolate.interp1d(
            zeta_t_rad, zeta_values, bounds_error=False, fill_value=np.nan
        )(t_rad)
        zeta = zeta.astype(float)

        if np.any(np.isnan(zeta)):
            warnings.warn(
                "t_rads outside of zeta factor interpolation"
                f" zeta_min={zeta_data.columns.values.min():.2f} zeta_max={zeta_data.columns.values.max():.2f} "
                "- replacing with 1s"
            )
            zeta[np.isnan(zeta)] = 1.0

        return zeta


class RadiationFieldCorrection(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    delta : Pandas DataFrame, dtype float
            Calculates the radiation field correction (see Mazzali & Lucy, 1993) if
            not given as input in the config. file. The default chi_0_species is
            Ca II, which is good for type Ia supernovae. For type II supernovae,
            (1, 1) should be used. Indexed by atomic number, ion number. The columns are zones.
    """

    outputs = ("delta",)
    latex_name = ("\\delta",)

    def __init__(
        self,
        plasma_parent=None,
        departure_coefficient=None,
        chi_0_species=(1, 1),
        delta_treatment=None,
    ):
        super(RadiationFieldCorrection, self).__init__(plasma_parent)
        self.departure_coefficient = departure_coefficient
        try:
            self.delta_treatment = self.plasma_parent.delta_treatment
        except:
            self.delta_treatment = delta_treatment

        self.chi_0_species = chi_0_species

    def _set_chi_0(self, ionization_data):
        if self.chi_0_species == (20, 2):
            self.chi_0 = 1.9020591570241798e-11
        else:
            self.chi_0 = ionization_data.loc[self.chi_0_species]

    def calculate(
        self, w, ionization_data, beta_rad, t_electrons, t_rad, beta_electron
    ):
        if getattr(self, "chi_0", None) is None:
            self._set_chi_0(ionization_data)
        if self.delta_treatment is None:
            if self.departure_coefficient is None:
                departure_coefficient = 1.0 / w
            else:
                departure_coefficient = self.departure_coefficient
            radiation_field_correction = -np.ones(
                (len(ionization_data), len(beta_rad))
            )
            less_than_chi_0 = (ionization_data < self.chi_0).values
            factor_a = t_electrons / (departure_coefficient * w * t_rad)
            radiation_field_correction[~less_than_chi_0] = factor_a * np.exp(
                np.outer(
                    ionization_data.values[~less_than_chi_0],
                    beta_rad - beta_electron,
                )
            )
            radiation_field_correction[less_than_chi_0] = 1 - np.exp(
                np.outer(
                    ionization_data.values[less_than_chi_0],
                    beta_rad,
                )
                - beta_rad * self.chi_0
            )
            radiation_field_correction[less_than_chi_0] += factor_a * np.exp(
                np.outer(
                    ionization_data.values[less_than_chi_0],
                    beta_rad,
                )
                - self.chi_0 * beta_electron
            )
        else:
            radiation_field_correction = (
                np.ones((len(ionization_data), len(beta_rad)))
                * self.plasma_parent.delta_treatment
            )
        delta = pd.DataFrame(
            radiation_field_correction,
            columns=np.arange(len(t_rad)),
            index=ionization_data.index,
        )
        return delta


class IonNumberDensity(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    ion_number_density : Pandas DataFrame, dtype float
                         Index atom number, ion number. Columns zones.
    electron_densities : Numpy Array, dtype float

    Convergence process to find the correct solution. A trial value for
    the electron density is initiated in a particular zone. The ion
    number densities are then calculated using the Saha equation. The
    electron density is then re-calculated by using the ion number
    densities to sum over the number of free electrons. If the two values
    for the electron densities are not similar to within the threshold
    value, a new guess for the value of the electron density is chosen
    and the process is repeated.
    """

    outputs = ("ion_number_density", "electron_densities")
    latex_name = (
        "N_{i,j}",
        "n_{e}",
    )

    def __init__(self, plasma_parent, ion_zero_threshold=1e-20):
        super(IonNumberDensity, self).__init__(plasma_parent)
        self.ion_zero_threshold = ion_zero_threshold
        self.block_ids = None

    def update_helium_nlte(self, ion_number_density, number_density):
        ion_number_density.loc[2].loc[0] = 0.0
        ion_number_density.loc[2].loc[2] = 0.0
        ion_number_density.loc[2].loc[1].update(number_density.loc[2])
        return ion_number_density

    def calculate_with_n_electron(
        self, phi, partition_function, number_density, n_electron
    ):
        if self.block_ids is None:
            self.block_ids = self._calculate_block_ids(phi)

        ion_populations = np.empty_like(partition_function.values)

        phi_electron = np.nan_to_num(phi.values / n_electron.values)

        for i, start_id in enumerate(self.block_ids[:-1]):
            end_id = self.block_ids[i + 1]
            current_phis = phi_electron[start_id:end_id]
            phis_product = np.cumprod(current_phis, 0)

            tmp_ion_populations = np.empty(
                (current_phis.shape[0] + 1, current_phis.shape[1])
            )
            tmp_ion_populations[0] = number_density.values[i] / (
                1 + np.sum(phis_product, axis=0)
            )
            tmp_ion_populations[1:] = tmp_ion_populations[0] * phis_product

            ion_populations[start_id + i : end_id + 1 + i] = tmp_ion_populations

        ion_populations[ion_populations < self.ion_zero_threshold] = 0.0

        return pd.DataFrame(
            data=ion_populations, index=partition_function.index
        )

    @staticmethod
    def _calculate_block_ids(phi):
        return calculate_block_ids_from_dataframe(phi)

    def calculate(self, phi, partition_function, number_density):
        print("WARNING, WARNING \n\n")
        print("WARNING, WARNING \n\n")
        print("WARNING, WARNING \n\n")
        n_e_convergence_threshold = 0.05
        n_electron = number_density.sum(axis=0)
        n_electron_iterations = 0

        while True:
            ion_number_density = self.calculate_with_n_electron(
                phi, partition_function, number_density, n_electron
            )
            if hasattr(self.plasma_parent, "plasma_properties_dict"):
                if (
                    "HeliumNLTE"
                    in self.plasma_parent.plasma_properties_dict.keys()
                ):
                    ion_number_density = self.update_helium_nlte(
                        ion_number_density, number_density
                    )
            ion_numbers = ion_number_density.index.get_level_values(1).values
            ion_numbers = ion_numbers.reshape((ion_numbers.shape[0], 1))
            new_n_electron = (ion_number_density.values * ion_numbers).sum(
                axis=0
            )
            if np.any(np.isnan(new_n_electron)):
                raise PlasmaIonizationError(
                    'n_electron just turned "nan" - aborting'
                )
            n_electron_iterations += 1
            if n_electron_iterations > 100:
                logger.warning(
                    f"n_electron iterations above 100 ({n_electron_iterations}) -"
                    " something is probably wrong"
                )
            if np.all(
                np.abs(new_n_electron - n_electron) / n_electron
                < n_e_convergence_threshold
            ):
                break
            n_electron = 0.5 * (new_n_electron + n_electron)

        return ion_number_density, n_electron


class LTEIonNumberDensity(IonNumberDensity):
    outputs = ("lte_ion_number_density",)
    latex_name = ("N_{i,j}^*",)

    def calculate(
        self,
        phi_Te,
        lte_partition_function_Te,
        number_density,
        electron_densities,
    ):
        return self.calculate_with_n_electron(
            phi_Te,
            lte_partition_function_Te,
            number_density,
            electron_densities,
        )


class NLTEIonNumberDensity(ProcessingPlasmaProperty):
    outputs = ("ion_number_density", "electron_densities")
    latex_name = (
        "N_{i,j}",
        "n_{e}",
    )

    def calculate(
        self,
        phi,
        alpha_stim,
        alpha_sp,
        gamma,
        coll_ion_coeff,
        coll_recomb_coeff,
        number_density,
        level_boltzmann_factor,
    ):
        partition_function = PartitionFunction.calculate(level_boltzmann_factor)
        ion_number_densities = pd.DataFrame(
            index=partition_function.index, columns=[]
        )
        electron_densities = np.zeros(partition_function.shape[1], dtype=float)
        alpha_tot = gamma_tot = coll_ion_tot = coll_recomb_tot = None

        number_density_vectors, number_density_vectors_neutral = (
            self._setup_number_density_vectors(
                number_density, partition_function.index
            )
        )
        charge_conservation_vector = self._get_charge_conservation_vector(
            partition_function.index
        )
        phi_prep = self._prepare_phi(phi, partition_function.index)
        number_conservation_index = self._get_number_conservation_index(
            partition_function.index
        )
        if gamma is not None:
            level_pop_fractions = self._calculate_level_pop_fractions(
                level_boltzmann_factor, partition_function
            )
            alpha_tot = self._calculate_alpha_tot(
                alpha_stim, alpha_sp, partition_function.index
            )
            gamma_tot = self._calculate_tot_ion_rate(
                gamma, level_pop_fractions, partition_function.index
            )

            coll_ion_tot = self._calculate_tot_ion_rate(
                coll_ion_coeff, level_pop_fractions, partition_function.index
            )
            coll_recomb_tot = self._calculate_coll_recomb_tot(
                coll_recomb_coeff, partition_function.index
            )
        for i in range(partition_function.shape[1]):
            initial_guesses = self._generate_initial_guesses(
                number_density_vectors,
                number_density_vectors_neutral,
                number_density,
                shell=i,
            )
            for initial_guess in initial_guesses:
                rates = [alpha_tot, gamma_tot, coll_ion_tot, coll_recomb_tot]
                func, jac = self._get_functions(
                    charge_conservation_vector,
                    number_density_vectors,
                    number_conservation_index,
                    phi_prep,
                    *rates,
                    shell=i,
                )
                sol = root(
                    fun=func, x0=initial_guess, jac=jac, options={"xtol": 1e-12}
                )
                ion_number_density = pd.Series(
                    sol.x[:-1], index=partition_function.index
                )

                ion_number_density = self._prepare_ion_densities(
                    ion_number_density
                )
                success_number_conservation = self._check_number_conservation(
                    ion_number_density, number_density[i]
                )
                success = sol.success and success_number_conservation
                if success:
                    break

            else:
                raise PlasmaIonizationError(
                    f"Ionization calculation failed in shell {i}."
                )

            ion_number_densities[i] = ion_number_density
            electron_densities[i] = sol.x[-1]

        if gamma is not None:
            # import pdb; pdb.set_trace()
            if self.plasma_parent.previous_ion_number_density is not None:
                conv = np.fabs(
                    (
                        self.plasma_parent.previous_ion_number_density
                        - ion_number_densities
                    )
                    / self.plasma_parent.previous_ion_number_density
                )
                # print "Ionization Convergence:", conv.loc[(1,0)]
                # print "Ionization Convergence max:", conv.max().max()
                # if conv.max().max() > 0.01:
                #    self.plasma_parent.update(previous_ion_number_density=ion_number_densities,
                #                          previous_electron_densities=pd.Series(electron_densities))

            # self.plasma_parent.previous_electron_densities.update(pd.Series(electron_densities))
            # self.plasma_parent.previous_ion_number_density.update(ion_number_densities)

        return ion_number_densities, pd.Series(electron_densities)

    @staticmethod
    def _prepare_ion_densities(ion_number_density):
        negative_elements = ion_number_density < 0.0
        ion_number_density[negative_elements] = 0.0

        ion_zero_threshold = 1e-20
        ion_number_density[ion_number_density < ion_zero_threshold] = 0.0
        return ion_number_density

    @staticmethod
    def _check_number_conservation(
        ion_number_density, number_density, rtol=1e-5
    ):
        number_density_from_ions = ion_number_density.groupby(level=0).sum()
        rdiff = np.fabs(number_density_from_ions - number_density).divide(
            number_density
        )
        success_number_conservation = (rdiff < rtol).all()
        if not success_number_conservation:
            logger.warning(f"Number conservation violated:\n{rdiff}")
        return success_number_conservation

    def _get_functions(
        self,
        charge_conservation_vector,
        number_density_vectors,
        number_conservation_index,
        phi_prep,
        alpha_tot,
        gamma_tot,
        coll_ion_tot,
        coll_recomb_tot,
        shell,
    ):
        lte_diag = -np.hstack([phi_prep[shell].values, np.zeros(1)])
        lte_offdiag = (lte_diag != 0).astype(float)[:-1]
        number_density = np.hstack([number_density_vectors[shell], np.zeros(1)])
        if alpha_tot is not None:
            alpha_matrix = self._setup_recomb_rate_matrix(
                alpha_tot[shell].values
            )
            gamma_matrix = self._setup_ion_rate_matrix(gamma_tot[shell].values)

            coll_recomb_matrix = self._setup_recomb_rate_matrix(
                coll_recomb_tot[shell].values
            )
            coll_ion_matrix = self._setup_ion_rate_matrix(
                coll_ion_tot[shell].values
            )

            nlte_mask = (
                alpha_matrix
                + gamma_matrix
                + coll_ion_matrix
                + coll_recomb_matrix
            ) != 0.0
            nlte_species_mask = self._get_nlte_species_mask(gamma_tot)
            nlte_mask = np.logical_and(nlte_mask, nlte_species_mask)

        def setup_rate_matrix(x):
            electron_density = x[-1]
            rate_matrix = (
                np.diag(lte_diag) + np.diag(lte_offdiag, k=1) * electron_density
            )
            if alpha_tot is not None:
                rate_matrix[nlte_mask] = 0.0
                rate_matrix[nlte_mask] += (
                    alpha_matrix * electron_density + gamma_matrix
                )[nlte_mask]
                rate_matrix[nlte_mask] += (
                    coll_ion_matrix * electron_density
                    + coll_recomb_matrix * electron_density**2.0
                )[nlte_mask]
            rate_matrix[number_conservation_index] = 1.0
            rate_matrix[-1] = charge_conservation_vector
            return rate_matrix

        def func(x):
            rate_matrix = setup_rate_matrix(x)
            return np.dot(rate_matrix, x) - number_density

        def jac(x):
            electron_density = x[-1]
            rate_matrix = setup_rate_matrix(x)

            derivative_matrix = np.diag(lte_offdiag, k=1)
            if alpha_tot is not None:
                derivative_matrix[nlte_mask] = 0.0
                derivative_matrix += alpha_matrix + coll_ion_matrix
                derivative_matrix += coll_recomb_matrix * electron_density
                derivative_matrix[number_conservation_index] = 0.0
            derivative_by_electron_density = np.dot(derivative_matrix, x)[:-1]
            rate_matrix[:-1, -1] = derivative_by_electron_density
            return rate_matrix

        return func, jac

    def _get_nlte_elements_mask(self, gamma_tot):
        nlte_elements = [
            species[0] for species in self.plasma_parent.nlte_species
        ]
        nlte_elements_mask = np.isin(
            gamma_tot.index.get_level_values(0), nlte_elements
        )
        nlte_elements_mask = np.hstack(
            [nlte_elements_mask, np.zeros(1, dtype=bool)]
        )
        nlte_elements_mask = nlte_elements_mask.reshape(1, -1)
        nlte_elements_mask = np.dot(nlte_elements_mask.T, nlte_elements_mask)
        return nlte_elements_mask

    def _get_nlte_species_mask(self, gamma_tot):
        nlte_species = self.plasma_parent.nlte_species

        nlte_species_mask = np.zeros(len(gamma_tot))
        for species in nlte_species:
            mask1 = [s == species for s in gamma_tot.index.values]
            nlte_species_mask = np.logical_or(nlte_species_mask, mask1)
        nlte_species_mask = np.hstack(
            [nlte_species_mask, np.zeros(1, dtype=bool)]
        )
        nlte_species_mask = nlte_species_mask.reshape(-1, 1)
        nlte_species_mask = np.tile(nlte_species_mask, len(nlte_species_mask))
        return nlte_species_mask

    @staticmethod
    def _setup_recomb_rate_matrix(recomb_rate):
        offdiag = recomb_rate
        diag = np.hstack([np.zeros(1), -recomb_rate])
        return np.diag(diag) + np.diag(offdiag, k=1)

    @staticmethod
    def _setup_ion_rate_matrix(ion_rate):
        offdiag = ion_rate
        diag = np.hstack([-ion_rate, np.zeros(1)])
        return np.diag(diag) + np.diag(offdiag, k=-1)

    @staticmethod
    def _prepare_phi(phi, ion_index):
        # Check for Nans
        no_nans = pd.isnull(phi).sum().sum()
        if no_nans:
            logger.warning(f"Phis contain {no_nans} Nans")
            phi = phi.fillna(phi.min().min())
        # Zero phi values pose a problem for the root finding algorithm. Set them to a small value.
        phi[phi == 0.0] = 1.0e-10 * phi[phi > 0.0].min().min()

        atomic_number = phi.index.get_level_values(0).values
        ion_number = phi.index.get_level_values(1).values
        new_index = pd.MultiIndex.from_arrays([atomic_number, ion_number - 1])
        phi_prep = phi.set_index(new_index).reindex(ion_index).fillna(0.0)
        return phi_prep

    @staticmethod
    def _get_number_conservation_index(ion_index):
        atomic_number = np.unique(
            ion_index.get_level_values(0).values.astype(int)
        )
        sum1 = (atomic_number + 1).cumsum() - 1
        index1 = np.concatenate(
            [np.ones(j + 1, dtype=int) * i for i, j in zip(sum1, atomic_number)]
        )
        index2 = np.arange(len(ion_index), dtype=int)
        return (index1, index2)

    @staticmethod
    def _setup_number_density_vectors(number_density, ion_index):
        atomic_number = number_density.index.values
        tmp_index = pd.MultiIndex.from_arrays([atomic_number, atomic_number])
        tmp_index_neutral = pd.MultiIndex.from_arrays(
            [atomic_number, np.zeros_like(atomic_number)]
        )

        number_density_vectors = (
            (number_density.set_index(tmp_index)).reindex(ion_index).fillna(0)
        )
        number_density_neutral_vectors = (
            (number_density.set_index(tmp_index_neutral))
            .reindex(ion_index)
            .fillna(0)
        )
        return number_density_vectors, number_density_neutral_vectors

    @staticmethod
    def _get_charge_conservation_vector(ion_index):
        ion_charge = ion_index.get_level_values(1)
        return np.hstack([ion_charge, np.array([-1.0])])

    def _generate_initial_guesses(
        self,
        number_conservation_vectors,
        number_density_neutral,
        number_density,
        shell,
    ):
        if self.plasma_parent.previous_ion_number_density is not None:
            initial_ion_density = (
                self.plasma_parent.previous_ion_number_density[shell].values
            )
            initial_electron_density = (
                self.plasma_parent.previous_electron_densities[shell]
            )
            yield np.hstack(
                [initial_ion_density, np.array(initial_electron_density)]
            )

        # First guess: Everything is singly ionized
        index = number_density_neutral.index
        initial_ion_density = pd.Series(
            np.zeros_like(number_conservation_vectors[shell].values),
            index=index,
        )
        singly_ionized_index = pd.MultiIndex.from_product(
            [number_density.index, [1]]
        ).values
        initial_ion_density.loc[singly_ionized_index] = number_density[
            shell
        ].values
        yield np.hstack([initial_ion_density, initial_ion_density.sum()])

        # Second guess: Everything is ionized.
        ion_number = number_conservation_vectors.index.get_level_values(
            1
        ).values
        initial_ion_density = number_conservation_vectors[shell].values
        initial_electron_density = (ion_number * initial_ion_density).sum()
        yield np.hstack(
            [initial_ion_density, np.array(initial_electron_density)]
        )

        # Third guess: Everything is neutral. A small amount of free electrons is added.
        logger.warning(
            "Calculation did not succeed with first guess (=everything ionized)"
        )
        neutral = number_density_neutral[shell].values
        yield np.hstack([neutral, np.ones(1) * neutral.sum() * 1.0e-5])

        # Fourth guess: Every ionization stage within an ion has the same number density.
        logger.warning(
            "Calculation did not succeed with second guess (=everything neutral)"
        )
        ion_guess = []
        for name, group in number_density_neutral[shell].groupby(level=[0]):
            count = len(group)
            total_density = group[0]
            ion_guess.append(np.ones(count) * total_density / count)
        ion_guess = np.hstack(ion_guess)
        ion_number = number_density_neutral.index.get_level_values(1).values
        electron_density = (ion_guess * ion_number).sum()
        yield np.hstack([ion_guess, np.array(electron_density)])

    # TODO: Probably only need this for nlte species
    def _calculate_level_pop_fractions(
        self, level_boltzmann_factor, partition_function
    ):
        boltzmann_factor = self._prepare_boltzmann_factor(
            level_boltzmann_factor
        )
        partition_function_index = get_ion_multi_index(
            boltzmann_factor.index, next_higher=False
        )
        partition_function = partition_function.loc[
            partition_function_index
        ].values
        return boltzmann_factor.divide(partition_function)

    @staticmethod
    def _prepare_boltzmann_factor(boltzmann_factor):
        atomic_number = boltzmann_factor.index.get_level_values(0)
        ion_number = boltzmann_factor.index.get_level_values(1)
        selected_ions_mask = atomic_number != ion_number
        return boltzmann_factor[selected_ions_mask]

    @staticmethod
    def _calculate_alpha_tot(alpha_stim, alpha_sp, index):
        alpha_tot = (alpha_sp + alpha_stim).groupby(level=[0, 1]).sum()
        return alpha_tot.reindex(index).fillna(0)

    @staticmethod
    def _calculate_coll_recomb_tot(coll_recomb_coeff, index):
        coll_recomb_tot = coll_recomb_coeff.groupby(level=[0, 1]).sum()
        return coll_recomb_tot.reindex(index).fillna(0)

    @staticmethod
    def _calculate_tot_ion_rate(rate, level_population_fractions, index):
        total_rate = (
            rate.multiply(level_population_fractions.loc[rate.index])
            .groupby(level=[0, 1])
            .sum()
        )
        return total_rate.reindex(index).fillna(0)

    # Not used atm. Intended as input for 'diag'-option of 'hybr' root finding method.
    @staticmethod
    def _get_scale_factors(number_density_vectors, shell):
        counts = number_density_vectors[0].groupby(level=0).count().values
        number_density = (
            number_density_vectors.groupby(level=0).max()[shell].values
        )
        scale_factor = []
        for i, count in enumerate(counts):
            scale_factor.append(np.ones(count, dtype=float) * number_density[i])
        scale_factor.append(np.array(number_density.sum()))
        scale_factor = np.hstack(scale_factor)
        scale_factor /= scale_factor.max()
        return 1.0 / scale_factor
