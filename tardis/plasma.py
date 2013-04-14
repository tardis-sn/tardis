#Calculations of the Plasma conditions

#import constants
import numpy as np
import logging
from astropy import table, units, constants

from pandas import DataFrame, Series, Index, lib as pdlib
import pandas as pd

import macro_atom

logger = logging.getLogger(__name__)
import pdb




#Defining soboleve constant

sobolev_coefficient = ((np.pi * constants.e.gauss.value ** 2) / (constants.m_e.cgs.value * constants.c.cgs.value))


class PlasmaException(Exception):
    pass


def intensity_black_body(nu, T):
    """
        Calculate the intensity of a black-body according to the following formula

        .. math::
            I(\\nu, T) = \frac{2h\\nu^3}{c^2}\frac{1}{e^{h\\nu \\beta_\\textrm{rad}} - 1}

    """
    beta_rad = 1 / (constants.k_B.cgs.value * T)

    return (2 * (constants.h.cgs.value * nu ** 3) / (constants.c.cgs.value ** 2)) / (
        np.exp(constants.h.cgs.value * nu * beta_rad) - 1)


class BasePlasma(object):
    """
    Model for BasePlasma

    Parameters
    ----------

    abundances : `~dict`
        A dictionary with the abundances for each element, e.g. {'Fe':0.5, 'Ni':0.5}

    t_rad : `~float`
        Temperature in Kelvin for the plasma

    density : `float`
        density in g/cm^3

        .. warning::
            Instead of g/cm^ will later use the keyword `density_unit` as unit

    atom_data : `~tardis.atomic.AtomData`-object

    max_ion_number : `~int`
        maximum used ionization of atom used in the calculation (inclusive the number)

    """

    @classmethod
    def from_abundance(cls, t_rad, w, abundance, density, atom_data, time_explosion, t_electron=None,
                       use_macro_atom=False, zone_id=None, nlte_species=[]):


    def __init__(self, t_rad, w, number_density, atom_data, time_explosion, t_electron=None, use_macro_atom=False,
                 zone_id=None, nlte_species=[]):
        self.number_density = number_density
        self.atom_data = atom_data
        self.initialize = True
        self.zone_id = zone_id
        self.time_explosion = time_explosion
        self.update_radiationfield(t_rad, w)
        self.nlte_species = nlte_species


    @property
    def t_rad(self):
        return self._t_rad

    @t_rad.setter
    def t_rad(self, value):
        self._t_rad = value
        self._beta_rad = 1 / (constants.k_B.cgs.value * self._t_rad)


    def validate_atom_data(self):
        required_attributes = ['lines', 'levels']
        for attribute in required_attributes:
            if not hasattr(self.atom_data, attribute):
                raise ValueError('AtomData incomplete missing')


    def update_radiationfield(self, t_rad, w):
        """
        This functions updates the radiation temperature `t_rad` and calculates the beta_rad
        Parameters
        ----------
        t_rad : float


        """
        self.t_rad = t_rad
        self.w = w


class LTEPlasma(BasePlasma):
    """
    Model for BasePlasma using a local thermodynamic equilibrium approximation.

    Parameters
    ----------

    abundances : `~dict`
       A dictionary with the abundances for each element

    t_rad : `~float`
       Temperature in Kelvin for the plasma

    density : `float`
       density in g/cm^3

       .. warning::
           Instead of g/cm^ will later use the keyword `density_unit` as unit

    atom_data : `~tardis.atomic.AtomData`-object

    """

    def __init__(self, number_density, atom_data, time_explosion, max_ion_number=None,
                 use_macro_atom=False, zone_id=None):
        BasePlasma.__init__(self, number_density, atom_data, time_explosion,
                            max_ion_number=max_ion_number, use_macro_atom=use_macro_atom, zone_id=zone_id)

        self.ion_number_density = None
        self.nlte_species = nlte_species


    def update_radiationfield(self, t_rad, t_electron=None, n_e_convergence_threshold=0.05, coronal_approximation=False,
                              classical_nebular=False):
        """
            This functions updates the radiation temperature `t_rad` and calculates the beta_rad
            Parameters. Then calculating :math:`g_e=\\left(\\frac{2 \\pi m_e k_\\textrm{B}T}{h^2}\\right)^{3/2}`.
            Next will calculate the partition functions, followed by the phis
            (using `calculate_saha`).

            Parameters
            ----------
            t_rad : float

            n_e_convergence_threshold : float
                The electron density convergence threshold. The number to stop when iterating over calculating the
                ionization balance.

       """
        BasePlasma.update_radiationfield(self, t_rad)

        if t_electron is None:
            self.t_electron = t_rad * 0.9
        else:
            self.t_electron = t_electron

        self.calculate_partition_functions(initialize=self.initialize)

        self.ge = ((2 * np.pi * constants.m_e.cgs.value / self.beta_rad) / (constants.h.cgs.value ** 2)) ** 1.5

        #Calculate the Saha ionization balance fractions
        phis = self.calculate_saha()

        #initialize electron density with the sum of number densities
        electron_density = self.number_density.sum()

        n_e_iterations = 0

        while True:
            self.calculate_ionization_balance(phis, electron_density)
            ion_numbers = np.array([item[1] for item in self.ion_number_density.index])
            new_electron_density = np.sum(self.ion_number_density.values * ion_numbers)
            n_e_iterations += 1
            if abs(new_electron_density - electron_density) / electron_density < n_e_convergence_threshold: break
            electron_density = 0.5 * (new_electron_density + electron_density)
        self.electron_density = new_electron_density
        logger.debug('Took %d iterations to converge on electron density' % n_e_iterations)

        self.calculate_level_populations()
        self.calculate_tau_sobolev()
        self.calculate_nlte_level_populations(coronal_approximation=coronal_approximation,
                                              classical_nebular=classical_nebular)

        if self.initialize:
            self.initialize = False

    def set_j_blues(self, j_blues=None):
        if j_blues is None:
            self.j_blues = intensity_black_body(self.atom_data.lines['nu'].values, self.t_rad)
        else:
            self.j_blues = j_blues


    def calculate_partition_functions(self, initialize=False):
        """
        Calculate partition functions for the ions using the following formula, where
        :math:`i` is the atomic_number, :math:`j` is the ion_number and :math:`k` is the level number.

        .. math::
            Z_{i,j} = \\sum_{k=0}^{max(k)_{i,j}} g_k \\times e^{-E_k / (k_\\textrm{b} T)}



        if self.initialize is True set the first time the partition functions are initialized.
        This will set a self.partition_functions and initialize with LTE conditions.


        Returns
        -------

        partition_functions : `~astropy.table.Table`
            with fields atomic_number, ion_number, partition_function

        """


        def group_calculate_partition_function(group):
            return np.sum(group['g'] *
                          np.exp(-group['energy'] * self.beta_rad))


        if self.initialize:
            logger.debug('Initializing the partition functions and indices')

            self.partition_functions = self.atom_data.levels.groupby(level=['atomic_number', 'ion_number']).apply(
                group_calculate_partition_function)

            self.atom_data.atom_ion_index = Series(np.arange(len(self.partition_functions)),
                                                   self.partition_functions.index)
            self.atom_data.levels_index2atom_ion_index = self.atom_data.atom_ion_index.ix[
                self.atom_data.levels.index.droplevel(2)].values
        else:
            if not hasattr(self, 'partition_functions'):
                raise ValueError("Called calculate partition_functions without initializing at least once")

            for species, group in self.atom_data.levels.groupby(level=['atomic_number', 'ion_number']):
                if species in self.nlte_species:
                    ground_level_population = self.level_populations[species][0]
                    self.partition_functions.ix[species] = self.atom_data.levels.ix[species]['g'][0] * \
                                                           np.sum(self.level_populations[
                                                                      species].values / ground_level_population)
                else:
                    self.partition_functions.ix[species] = np.sum(group['g'] * np.exp(-group['energy'] * self.beta_rad))


    def calculate_saha(self):
        """
        Calculating the ionization equilibrium using the Saha equation, where i is atomic number,
        j is the ion_number, :math:`n_e` is the electron density, :math:`Z_{i, j}` are the partition functions
        and :math:`\chi` is the ionization energy.

        .. math::


            \\Phi_{i,j} = \\frac{N_{i, j+1} n_e}{N_{i, j}}

            \\Phi_{i, j} = g_e \\times \\frac{Z_{i, j+1}}{Z_{i, j}} e^{-\chi_{j\\rightarrow j+1}/k_\\textrm{B}T}

        """

        def calculate_phis(group):
            return group[1:] / group[:-1].values

        phis = self.partition_functions.groupby(level='atomic_number').apply(calculate_phis)

        phis = Series(phis.values, phis.index.droplevel(0))

        phis *= self.ge * np.exp(-self.beta_rad * self.atom_data.ionization_data.ix[phis.index]['ionization_energy'])

        return phis


    def calculate_ionization_balance(self, phis, electron_density):
        """
        Calculate the ionization balance

        .. math::
            N(X) = N_1 + N_2 + N_3 + \\dots

            N(X) = (N_2/N_1) \\times N_1 + (N3/N2) \\times (N_2/N_1) \\times N_1 + \\dots

            N(X) = N_1(1 + N_2/N_1 + (N_3/N_2) \\times (N_2/N_1) + \\dots

            N(X) = N_1(1+ \\Phi_{i,j}/N_e + \\Phi_{i, j}/N_e \\times \\Phi_{i, j+1}/N_e + \\dots)


        """

        if self.ion_number_density is None:
            self.ion_number_density = pd.Series(index=self.partition_functions.index.copy())
            self.cleaned_levels = pd.Series(index=self.partition_functions.index.copy())

        for atomic_number, groups in phis.groupby(level='atomic_number'):
            current_phis = groups.values / electron_density
            phis_product = np.cumproduct(current_phis)

            neutral_atom_density = self.number_density.ix[atomic_number] / (1 + np.sum(phis_product))
            ion_densities = [neutral_atom_density] + list(neutral_atom_density * phis_product)

            self.ion_number_density.ix[atomic_number] = ion_densities


    def calculate_level_populations(self):
        """
        Calculate the level populations and storing in self.level_populations table.
        :math:`N` denotes the ion number density calculated with `calculate_ionization_balance`, i is the atomic number,
        j is the ion number and k is the level number.

        .. math::
            \\frac{g_k}{Z_{i,j}} \\times N_{i, j} \\times e^{-\\beta_\\textrm{rad} \\times E_k}

        """

        Z = self.partition_functions.values[self.atom_data.levels_index2atom_ion_index]

        ion_number_density = self.ion_number_density.values[self.atom_data.levels_index2atom_ion_index]

        levels_g = self.atom_data.levels['g'].values
        levels_energy = self.atom_data.levels['energy'].values
        level_populations = (levels_g / Z) * ion_number_density * np.exp(-self.beta_rad * levels_energy)

        if self.initialize:
            self.level_populations = Series(level_populations, index=self.atom_data.levels.index)

        else:
            level_populations = Series(level_populations, index=self.atom_data.levels.index)
            self.level_populations.update(level_populations[~self.atom_data.nlte_data.nlte_mask])


    def calculate_nlte_level_populations(self, coronal_approximation=False, classical_nebular=False):
        """
        Calculating the NLTE level populations for specific ions

        """

        if not hasattr(self, 'beta_sobolevs'):
            self.beta_sobolevs = np.zeros_like(self.atom_data.lines['nu'].values)

        macro_atom.calculate_beta_sobolev(self.tau_sobolevs, self.beta_sobolevs)

        if coronal_approximation:
            beta_sobolevs = np.ones_like(self.beta_sobolevs)
            j_blues = np.zeros_like(self.j_blues)
        else:
            beta_sobolevs = self.beta_sobolevs
            j_blues = self.j_blues

        if classical_nebular:
            print "setting classical nebular = True"
            beta_sobolevs[:] = 1.0

        for species in self.nlte_species:
            logger.info('Calculating rates for species %s', species)
            number_of_levels = self.level_populations.ix[species].size

            level_populations = self.level_populations.ix[species].values
            lnl = self.atom_data.nlte_data.lines_level_number_lower[species]
            lnu = self.atom_data.nlte_data.lines_level_number_upper[species]

            lines_index = self.atom_data.nlte_data.lines_idx[species]
            A_uls = self.atom_data.nlte_data.A_uls[species]
            B_uls = self.atom_data.nlte_data.B_uls[species]
            B_lus = self.atom_data.nlte_data.B_lus[species]

            r_lu_index = lnu * number_of_levels + lnl
            r_ul_index = lnl * number_of_levels + lnu

            r_ul_matrix = np.zeros((number_of_levels, number_of_levels), dtype=np.float64)
            r_ul_matrix.ravel()[r_ul_index] = A_uls
            r_ul_matrix.ravel()[r_ul_index] *= beta_sobolevs[lines_index]

            stimulated_emission_matrix = np.zeros_like(r_ul_matrix)
            stimulated_emission_matrix.ravel()[r_lu_index] = 1 - (level_populations[lnu] * B_uls) / (
                level_populations[lnl] * B_lus)

            #stimulated_emission_matrix[stimulated_emission_matrix < 0.] = 0.0

            r_lu_matrix = np.zeros_like(r_ul_matrix)
            r_lu_matrix.ravel()[r_lu_index] = B_lus * j_blues[lines_index] * beta_sobolevs[lines_index]
            r_lu_matrix *= stimulated_emission_matrix

            collision_matrix = self.atom_data.nlte_data.get_collision_matrix(species,
                                                                             self.t_electron) * self.electron_density

            rates_matrix = r_lu_matrix + r_ul_matrix + collision_matrix

            for i in xrange(number_of_levels):
                rates_matrix[i, i] = -np.sum(rates_matrix[:, i])

            rates_matrix[0] = 1.0

            x = np.zeros(rates_matrix.shape[0])
            x[0] = 1.0
            relative_level_populations = np.linalg.solve(rates_matrix, x)

            self.level_populations.ix[species] = relative_level_populations * self.ion_number_density.ix[species]

            return


"""
            #Cleaning Level populations
            self.cleaned_levels.ix[species] = 0
            self.lowest_cleaned_level = 100000
            for i in xrange(1, number_of_levels):
                n_upper = self.level_populations.ix[species][i]
                n_lower = self.level_populations.ix[species][i - 1]

                g_upper = float(self.atom_data.levels.ix[species]['g'][i])
                g_lower = float(self.atom_data.levels.ix[species]['g'][i - 1])

                current_stim_ems = (n_upper / n_lower) * (g_lower / g_upper)

                if current_stim_ems > 1.:
                    self.level_populations.ix[species[0], species[1], i] = (1 - 1e-12) * (g_upper / g_lower) * n_lower
                    self.cleaned_levels.ix[species] += 1
                    self.lowest_cleaned_level = min(i, self.lowest_cleaned_level)
                    #### After cleaning check if the normalization is good:

            if abs((self.level_populations.ix[species].sum() / self.ion_number_density.ix[species] - 1)) > 0.02:
                logger.warn("NLTE populations (after cleaning) does not sum up to 1 within 2 percent "
                            "(%.2f / 1.0 - zone id = %s, lowest_cleaned_level=%d)",
                            ((1 - self.level_populations.ix[species].sum() / self.ion_number_density.ix[species])),
                            self.zone_id, self.lowest_cleaned_level)

            logger.debug('Number of cleaned levels %d of %d (zone id =%s)', self.cleaned_levels.ix[species],
                         self.level_populations.ix[species].count(), self.zone_id)

            if float(self.cleaned_levels.ix[species]) / self.level_populations.ix[species].count() > 0.5:
                logger.warn('Number of cleaned levels very high %d of %d (zone id=%s, lowest_cleaned_level=%d)',
                            self.cleaned_levels.ix[species],
                            self.level_populations.ix[species].count(), self.zone_id, self.lowest_cleaned_level)
"""


def calculate_tau_sobolev(self):
    """
    This function calculates the Sobolev optical depth :math:`\\tau_\\textrm{Sobolev}`



    .. math::
        C_\\textrm{Sobolev} = \\frac{\\pi e^2}{m_e c}

        \\tau_\\textrm{Sobolev} = C_\\textrm{Sobolev}\,  \\lambda\\, f_{\\textrm{lower}\\rightarrow\\textrm{upper}}\\,
            t_\\textrm{explosion}\, N_\\textrm{lower}



    .. note::
        Currently we're ignoring the term for stimulated emission:
            :math:`(1 - \\frac{g_\\textrm{lower}}{g_\\textrm{upper}}\\frac{N_\\textrm{upper}}{N_\\textrm{lower}})`


    """

    f_lu = self.atom_data.lines['f_lu'].values
    wavelength = self.atom_data.lines['wavelength_cm'].values
    n_lower = self.level_populations.values[self.atom_data.lines_lower2level_idx]
    self.tau_sobolevs = sobolev_coefficient * f_lu * wavelength * self.time_explosion * n_lower


def calculate_bound_free(self):
    """
    :return:
    """
    nu_bins = range(1000, 10000, 1000) #TODO: get the binning from the input file.
    try:
        bf = np.zeros(len(self.atom_data.levels), len(self.atom_data.selected_atomic_numbers), len(nu_bins))
    except AttributeError:
        logger.critical("Err creating the bf array.")

    phis = self.calculate_saha()
    nnlevel = self.level_populations
    for nu in nu_bins:
        for i, (level_id, level) in enumerate(self.atom_data.levels.iterrows()):
            atomic_number = level.name[0]
            ion_number = level.name[1]
            level_number = level.name[2]
            sigma_bf_th = self.atom_data.ion_cx_th.ix[atomic_number, ion_number, level_number]
            phi = phis.ix[atomic_number, ion_number]


def update_macro_atom(self):
    """
        Updating the Macro Atom computations

    """

    macro_tau_sobolevs = self.tau_sobolevs[self.atom_data.macro_atom_data['lines_idx'].values.astype(int)]

    beta_sobolevs = np.zeros_like(macro_tau_sobolevs)

    macro_atom.calculate_beta_sobolev(macro_tau_sobolevs, beta_sobolevs)

    transition_probabilities = self.atom_data.macro_atom_data['transition_probability'] * beta_sobolevs

    transition_up_filter = self.atom_data.macro_atom_data['transition_type'] == 1

    j_blues = self.j_blues[self.atom_data.macro_atom_data['lines_idx'].values[transition_up_filter.__array__()]]

    transition_probabilities[transition_up_filter.__array__()] *= j_blues

    reference_levels = np.hstack((0, self.atom_data.macro_atom_references['count_total'].__array__().cumsum()))

    #Normalizing the probabilities
    #TODO speedup possibility save the new blockreferences with 0 and last block
    block_references = np.hstack((self.atom_data.macro_atom_references['block_references'],
                                  len(self.atom_data.macro_atom_data)))
    macro_atom.normalize_transition_probabilities(transition_probabilities, block_references)

    return transition_probabilities


class NebularPlasma(LTEPlasma):
    """
    Model for BasePlasma using the Nebular approximation

    Parameters
    ----------

    abundances : `~dict`
       A dictionary with the abundances for each element

    t_rad : `~float`
       Temperature in Kelvin for the plasma

    density : `float`
       density in g/cm^3

       .. warning::
           Instead of g/cm^ will later use the keyword `density_unit` as unit

    atom_data : `~tardis.atomic.AtomData`-object

    t_electron : `~float`, or `None`
        the electron temperature. if set to `None` we assume the electron temperature is 0.9 * radiation temperature

    """

    def __init__(self, number_density, atom_data, time_explosion, nlte_species=[], t_electron=None,
                 density_unit='g/cm^3',
                 max_ion_number=None,
                 use_macro_atom=False, zone_id=None):
        BasePlasma.__init__(self, number_density, atom_data, time_explosion=time_explosion, density_unit=density_unit,
                            max_ion_number=max_ion_number,
                            use_macro_atom=use_macro_atom, zone_id=zone_id)

        self.ion_number_density = None
        self.nlte_species = nlte_species


    def update_radiationfield(self, t_rad, w, t_electron=None, n_e_convergence_threshold=0.05,
                              coronal_approximation=False, classical_nebular=True):
        BasePlasma.update_radiationfield(self, t_rad)

        self.w = w

        if t_electron is None:
            self.t_electron = 0.9 * self.t_rad

        self.beta_electron = 1 / (self.t_electron * constants.k_B.cgs.value)

        self.calculate_partition_functions()

        self.ge = ((2 * np.pi * constants.m_e.cgs.value / self.beta_rad) / (constants.h.cgs.value ** 2)) ** 1.5
        #Calculate the Saha ionization balance fractions
        phis = self.calculate_saha()

        #initialize electron density with the sum of number densities
        electron_density = self.number_density.sum()

        n_e_iterations = 0

        while True:
            self.calculate_ionization_balance(phis, electron_density)
            ion_numbers = np.array([item[1] for item in self.ion_number_density.index])
            new_electron_density = np.sum(self.ion_number_density.values * ion_numbers)
            n_e_iterations += 1
            if abs(new_electron_density - electron_density) / electron_density < n_e_convergence_threshold: break
            electron_density = 0.5 * (new_electron_density + electron_density)

        self.electron_density = new_electron_density
        logger.debug('Took %d iterations to converge on electron density' % n_e_iterations)

        self.calculate_level_populations()
        self.calculate_tau_sobolev()
        self.calculate_nlte_level_populations(coronal_approximation=coronal_approximation,
                                              classical_nebular=classical_nebular)

        if self.initialize:
            self.initialize = False

    def set_j_blues(self, j_blues=None):
        if j_blues is None:
            self.j_blues = self.w * intensity_black_body(self.atom_data.lines['nu'].values, self.t_rad)
        else:
            self.j_blues = j_blues

    def calculate_partition_functions(self):
        """
        Calculate partition functions for the ions using the following formula, where
        :math:`i` is the atomic_number, :math:`j` is the ion_number and :math:`k` is the level number.

        .. math::
            Z_{i,j} = \\underbrace{\\sum_{k=0}^{max(k)_{i,j}} g_k \\times e^{-E_k / (k_\\textrm{b} T)}}_\\textrm{metastable levels} +
                    \\underbrace{W\\times\\sum_{k=0}^{max(k)_{i,j}} g_k \\times e^{-E_k / (k_\\textrm{b} T)}}_\\textrm{non-metastable levels}



        Returns
        -------

        partition_functions : `~astropy.table.Table`
            with fields atomic_number, ion_number, partition_function

        """

        def group_calculate_partition_function(group):
            metastable = group['metastable']
            meta_z = np.sum(group['g'][metastable] * np.exp(-group['energy'][metastable] * self.beta_rad))
            non_meta_z = np.sum(group['g'][~metastable] * np.exp(-group['energy'][~metastable] * self.beta_rad))
            return meta_z + self.w * non_meta_z


        if self.initialize:
            logger.debug('Initializing the partition functions and indices')

            self.partition_functions = self.atom_data.levels.groupby(level=['atomic_number', 'ion_number']).apply(
                group_calculate_partition_function)

            self.atom_data.atom_ion_index = Series(np.arange(len(self.partition_functions)),
                                                   self.partition_functions.index)
            self.atom_data.levels_index2atom_ion_index = self.atom_data.atom_ion_index.ix[
                self.atom_data.levels.index.droplevel(2)].values
        else:
            if not hasattr(self, 'partition_functions'):
                raise ValueError("Called calculate partition_functions without initializing at least once")

            for species, group in self.atom_data.levels.groupby(level=['atomic_number', 'ion_number']):
                if species in self.nlte_species:
                    ground_level_population = self.level_populations[species][0]
                    self.partition_functions.ix[species] = self.atom_data.levels.ix[species]['g'][0] * \
                                                           np.sum(self.level_populations[
                                                                      species].values / ground_level_population)
                else:
                    self.partition_functions.ix[species] = np.sum(group['g'] * np.exp(-group['energy'] * self.beta_rad))


    def calculate_saha(self):
        """
        Calculating the ionization equilibrium using the Saha equation, where i is atomic number,
        j is the ion_number, :math:`n_e` is the electron density, :math:`Z_{i, j}` are the partition functions
        and :math:`\chi` is the ionization energy. For the `NebularPlasma` we first calculate the
        ionization balance assuming LTE conditions (:math:`\\Phi_{i, j}(\\textrm{LTE})`) and use factors to more accurately
        describe the plasma. The two important factors are :math:`\\zeta` - a correction factor to take into account
        ionizations from excited states. The second factor is :math:`\\delta` , adjusting the ionization balance for the fact that
        there's more line blanketing in the blue.

        The :math:`\\zeta` factor for different temperatures is read in to the `~tardis.atomic.NebularAtomData` and then
        interpolated for the current temperature.

        The :math:`\\delta` factor is calculated with :method:`calculate_radiation_field_correction`.

        Finally the ionization balance is adjusted (as equation 14 in :cite:`1993A&A...279..447M`):

        .. math::


            \\Phi_{i,j} =& \\frac{N_{i, j+1} n_e}{N_{i, j}} \\\\

            \\Phi_{i, j} =& W \\times[\\delta \\zeta + W ( 1 - \\zeta)] \\left(\\frac{T_\\textrm{e}}{T_\\textrm{R}}\\right)^{1/2}
            \\Phi_{i, j}(\\textrm{LTE})

        """

        phis = super(NebularPlasma, self).calculate_saha()

        delta = self.calculate_radiation_field_correction()

        zeta = Series(index=phis.index)

        for idx in zeta.index:
            try:
                current_zeta = self.atom_data.zeta_data[idx](self.t_rad)
            except KeyError:
                current_zeta = 1.0

            zeta.ix[idx] = current_zeta

        phis *= self.w * (delta.ix[phis.index] * zeta + self.w * (1 - zeta)) * \
                (self.t_electron / self.t_rad) ** .5

        return phis


    def calculate_radiation_field_correction(self, departure_coefficient=None,
                                             chi_threshold_species=(20, 1)):
        """
        Calculating radiation field correction factors according to Mazzali & Lucy 1993 (:cite:`1993A&A...279..447M`; henceforth ML93)


        In ML93 the radiation field correction factor is denoted as :math:`\\delta` and is calculated in Formula 15 & 20

        The radiation correction factor changes according to a ionization energy threshold :math:`\\chi_\\textrm{T}`
        and the species ionization threshold (from the ground state) :math:`\\chi_0`.

        For :math:`\\chi_\\textrm{T} \\ge \\chi_0`

        .. math::
            \\delta = \\frac{T_\\textrm{e}}{b_1 W T_\\textrm{R}} \\exp(\\frac{\\chi_\\textrm{T}}{k T_\\textrm{R}} -
            \\frac{\\chi_0}{k T_\\textrm{e}})

        For :math:`\\chi_\\textrm{T} < \\chi_0`

        .. math::
            \\delta = 1 - \\exp(\\frac{\\chi_\\textrm{T}}{k T_\\textrm{R}} - \\frac{\\chi_0}{k T_\\textrm{R}}) + \\frac{T_\\textrm{e}}{b_1 W T_\\textrm{R}} \\exp(\\frac{\\chi_\\textrm{T}}{k T_\\textrm{R}} -
            \\frac{\\chi_0}{k T_\\textrm{e}}),

        where :math:`T_\\textrm{R}` is the radiation field Temperature, :math:`T_\\textrm{e}` is the electron temperature and W is the
        dilution factor.

        Parameters
        ----------
        phi_table : `~astropy.table.Table`
            a table containing the field 'atomic_number', 'ion_number', 'phi'

        departure_coefficient : `~float` or `~None`, optional
            departure coefficient (:math:`b_1` in ML93) For the default (`None`) it is set to 1/W.

        chi_threshold_species : `~tuple`, optional
            This describes which ionization energy to use for the threshold. Default is Calcium II
            (1044 Angstrom; useful for Type Ia)
            For Type II supernovae use Lyman break (912 Angstrom) or (1,1) as the tuple

        Returns
        -------

        This function adds a field 'delta' to the phi table given to the function

        """
        #factor delta ML 1993
        if departure_coefficient is None:
            departure_coefficient = 1 / float(self.w)

        chi_threshold = self.atom_data.ionization_data['ionization_energy'].ix[chi_threshold_species]

        radiation_field_correction = (self.t_electron / (departure_coefficient * self.w * self.t_rad)) * \
                                     np.exp(self.beta_rad * chi_threshold - self.beta_electron *
                                            self.atom_data.ionization_data['ionization_energy'])

        less_than_chi_threshold = self.atom_data.ionization_data['ionization_energy'] < chi_threshold

        radiation_field_correction[less_than_chi_threshold] += 1 - \
                                                               np.exp(self.beta_rad * chi_threshold - self.beta_rad *
                                                                      self.atom_data.ionization_data[
                                                                          less_than_chi_threshold]['ionization_energy'])

        return radiation_field_correction


    def calculate_level_populations(self):
        """
        Calculate the level populations and putting them in the column 'number-density' of the self.levels table.
        :math:`N` denotes the ion number density calculated with `calculate_ionization_balance`, i is the atomic number,
        j is the ion number and k is the level number. For non-metastable levels we add the dilution factor (W) to the calculation.

        .. math::

            N_{i, j, k}(\\textrm{metastable}) &= \\frac{g_k}{Z_{i, j}}\\times N_{i, j} \\times e^{-\\beta_\\textrm{rad} E_k} \\\\
            N_{i, j, k}(\\textrm{not metastable}) &= W\\frac{g_k}{Z_{i, j}}\\times N_{i, j} \\times e^{-\\beta_\\textrm{rad} E_k} \\\\


        This function updates the 'number_density' column on the levels table (or adds it if non-existing)
        """
        Z = self.partition_functions.values[self.atom_data.levels_index2atom_ion_index]

        ion_number_density = self.ion_number_density.values[self.atom_data.levels_index2atom_ion_index]

        levels_g = self.atom_data.levels['g'].values
        levels_energy = self.atom_data.levels['energy'].values
        level_populations = (levels_g / Z) * ion_number_density * np.exp(-self.beta_rad * levels_energy)


        #only change between lte plasma and nebular
        level_populations[~self.atom_data.levels['metastable']] *= self.w

        if self.initialize:
            self.level_populations = Series(level_populations, index=self.atom_data.levels.index)

        else:
            level_populations = Series(level_populations, index=self.atom_data.levels.index)
            self.level_populations.update(level_populations[~self.atom_data.nlte_data.nlte_mask])


