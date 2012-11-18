#Calculations of the Plasma conditions

#import constants
import numpy as np
import logging
from astropy import table, units, constants
from collections import OrderedDict

logger = logging.getLogger(__name__)

#Bnu = lambda nu, t: (2 * constants.h * nu ** 3 / constants.c ** 2) * np.exp(
#    1 / ((constants.h * nu) / (constants.kb * t)))

class Plasma(object):
    """
    Model for Plasma

    Parameters
    ----------

    abundances : `dict`
        A dictionary with the abundances for each element

    t_rad : `~float`
        Temperature in Kelvin for the plasma

    density : `float`
        density in g/cm^3

        .. warning::
            Instead of g/cm^ will later use the keyword `density_unit` as unit

    atom_data : `~tardis.atomic.AtomData`-object

    """

    #TODO make density a astropy.quantity
    def __init__(self, abundances, t_rad, density, atom_data, density_unit='g/cm^3'):
        self.atom_data = atom_data
        self.density = density

        self.abundances = self.calculate_atom_number_densities(abundances)

        #Filtering the levels & lines data
        levels_filter = np.zeros(len(self.atom_data._levels)).astype(np.bool)
        lines_filter = np.zeros(len(self.atom_data._lines)).astype(np.bool)

        for atomic_number in self.abundances['atomic_number']:
            levels_filter = levels_filter | (self.atom_data._levels['atomic_number'] == atomic_number)
            lines_filter = lines_filter | (self.atom_data._lines['atomic_number'] == atomic_number)

        self.levels = self.atom_data._levels[levels_filter]
        self.lines = self.atom_data._lines[lines_filter]

        self.lines['wavelength'].convert_units_to('cm')

        level_tuple = [tuple(item) for item in self.levels.__array__()[['atomic_number', 'ion_number', 'level_number']]]
        self.levels_dict = dict(zip(level_tuple, np.arange(len(self.levels))))
        self.update_t_rad(t_rad)


    def calculate_atom_number_densities(self, abundances):
        """
        Calculates the atom number density, using the following formula, where Z is the atomic number
        and X is the abundance fraction

        .. math::
            N_{Z} = \\frac{\\rho_\\textrm{total}\\times \\textrm{X}_\\textrm{Z}}{m_\\textrm{Z}}

        """

        #Converting abundances
        abundance_table = zip(
            *[(self.atom_data.symbol2atomic_number[key], value) for key, value in abundances.items()])
        abundance_table = table.Table(abundance_table, names=('atomic_number', 'abundance_fraction'))


        #Normalizing Abundances

        abundance_sum = abundance_table['abundance_fraction'].sum()

        if abs(abundance_sum - 1) > 1e-5:
            logger.warn('Abundances do not add up to 1 (Sum = %.4f). Renormalizing', (abundance_sum))

        abundance_table['abundance_fraction'] /= abundance_sum

        atom_masses = np.array(
            [self.atom_data._atom['mass'][atomic_number - 1] for atomic_number in abundance_table['atomic_number']])
        number_density = (self.density * abundance_table['abundance_fraction']) / atom_masses

        number_density_col = table.Column('number_density', number_density, units=units.Unit(1))
        abundance_table.add_column(number_density_col)

        return abundance_table

    def update_t_rad(self, t_rad):
        """
        This functions updates the radiation temperature `t_rad` and calculates the beta_rad
        Parameters
        ----------
        t_rad : float


        """
        self.t_rad = t_rad
        self.beta_rad = 1 / (constants.cgs.k_B.real * t_rad)


class LTEPlasma(Plasma):
    def __init__(self, abundances, t_rad, density, atom_data, density_unit='g/cm^3'):
        Plasma.__init__(self, abundances, t_rad, density, atom_data, density_unit=density_unit)

        self.update_t_rad(t_rad)


    def update_t_rad(self, t_rad, n_e_convergence_threshold=0.05):
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
        Plasma.update_t_rad(self, t_rad)

        self.partition_functions = self.calculate_partition_functions()

        self.ge = ((2 * np.pi * constants.cgs.m_e / self.beta_rad) / (constants.cgs.h ** 2)) ** 1.5

        #Calculate the Saha ionization balance fractions
        phis = self.calculate_saha()

        #initialize electron density with the sum of number densities
        electron_density = np.sum(self.abundances['number_density'])

        n_e_iterations = 0

        while True:
            ionization_balance, new_electron_density = self.calculate_ionization_balance(phis, electron_density)
            n_e_iterations += 1
            if abs(new_electron_density - electron_density) / electron_density < n_e_convergence_threshold: break
            electron_density = 0.5 * (new_electron_density + electron_density)

        logger.info('Took %d iterations to converge on electron density' % n_e_iterations)

        self.ion_densities = ionization_balance
        self.calculate_level_populations()


    def calculate_partition_functions(self):
        """
        Calculate partition functions for the ions using the following formula, where
        :math:`i` is the atomic_number, :math:`j` is the ion_number and :math:`k` is the level number.

        .. math::
            Z_{i,j} = \\sum_{k=0}^{max(k)_{i,j}} g_k \\times e^{-E_k / (k_\\textrm{b} T)}



        Returns
        -------

        partition_functions : `~astropy.table.Table`
            with fields atomic_number, ion_number, partition_function

        """

        unique_atom_ion = np.unique(self.levels.__array__()[['atomic_number', 'ion_number']])

        partition_table = []
        for atomic_number, ion_number in unique_atom_ion:
            levels = self.atom_data.get_levels(atomic_number, ion_number)
            partition_function = np.sum(levels['g'] * np.exp(-levels['energy'] * self.beta_rad))
            partition_table.append((atomic_number, ion_number, partition_function))

        partition_table = np.array(partition_table, dtype=[('atomic_number', np.int),
                                                           ('ion_number', np.int),
                                                           ('partition_function', np.float)])

        partition_table = table.Table(partition_table)
        return partition_table

    def calculate_saha(self):
        """
        Calculating the ionization equilibrium using the Saha equation, where i is atomic number,
        j is the ion_number, :math:`n_e` is the electron density, :math:`Z_{i, j}` are the partition functions
        and :math:`\chi` is the ionization energy.

        .. math::


            \\Phi_{i,j} = \\frac{N_{i, j+1} n_e}{N_{i, j}}

            \\Phi_{i, j} = g_e \\times \\frac{Z_{i, j+1}}{Z_{i, j}} e^{-\chi_{j\\rightarrow j+1}/k_\\textrm{B}T}

        """
        denominators = []
        numerators = []
        phis = []
        atomic_numbers = []
        for atomic_number in self.abundances['atomic_number']:
            partition_functions_filter = self.partition_functions['atomic_number'] == atomic_number
            current_partition_functions = self.partition_functions[partition_functions_filter]

            phi_value = self.ge * current_partition_functions['partition_function'][1:] /\
                        current_partition_functions['partition_function'][:-1]

            phi_value *= np.exp(-self.beta_rad *\
                                self.atom_data.get_ions(atomic_number)['ionization_energy'][:len(phi_value)])
            phis += list(phi_value)
            denominators += list(current_partition_functions['ion_number'][:-1])
            numerators += list(current_partition_functions['ion_number'][1:])
            atomic_numbers += [atomic_number] * len(phi_value)



        #phi = (n_j+1 * ne / nj)
        table_columns = []
        table_columns.append(table.Column('atomic_number', atomic_numbers))
        table_columns.append(table.Column('numerator', numerators))
        table_columns.append(table.Column('denominator', denominators))
        table_columns.append(table.Column('phi', phis))

        phi_table = table.Table(table_columns)

        return phi_table


    def calculate_ionization_balance(self, phis, electron_density):
        """
        Calculate the ionization balance

        .. math::
            N(X) = N_1 + N_2 + N_3 + \\dots

            N(X) = (N_2/N_1) \\times N_1 + (N3/N2) \\times (N_2/N_1) \\times N_1 + \\dots

            N(X) = N1(1 + N_2/N_1 + (N_3/N_2) \\times (N_2/N_1) + \\dots


        """

        atomic_numbers = []
        ion_numbers = []
        ion_densities = []

        new_electron_density = 0

        for atomic_number, abundance_fraction, number_density in self.abundances:
            atomic_number_filter = phis['atomic_number'] == atomic_number
            current_phis = phis[atomic_number_filter]
            phis_product = np.cumprod(current_phis['phi'] / electron_density)
            neutral_atom_density = number_density / (1 + np.sum(phis_product))

            atomic_numbers += [atomic_number] * (len(current_phis) + 1)
            ion_numbers += [0] + list(current_phis['numerator'])

            ion_densities += [neutral_atom_density] + list(neutral_atom_density * phis_product)

        table_columns = []
        table_columns.append(table.Column('atomic_number', atomic_numbers))
        table_columns.append(table.Column('ion_number', ion_numbers))
        table_columns.append(table.Column('number_density', ion_densities))

        ion_densities = table.Table(table_columns)

        new_electron_density = np.sum(ion_densities['number_density'] * ion_densities['ion_number'])

        return ion_densities, new_electron_density


    def calculate_level_populations(self):
        if 'number_density' not in self.levels.dtype.names:
            n_levels = table.Column('number_density', np.empty_like(self.levels['atomic_number']).astype(np.float),
                dtype=np.float)
            self.levels.add_column(n_levels)

        n_levels = self.levels['number_density']

        for atomic_number, ion_number, partition_function in self.partition_functions:
            ion_density_filter = self.ion_densities['atomic_number'] == atomic_number
            ion_density_filter = ion_density_filter & (self.ion_densities['ion_number'] == ion_number)

            current_ion_density = self.ion_densities[ion_density_filter]['number_density'][0]

            level_filter = self.levels['atomic_number'] == atomic_number
            level_filter = level_filter & (self.levels['ion_number'] == ion_number)

            current_levels = self.levels[level_filter]

            n_levels[level_filter] = (current_levels['g'] / partition_function) * current_ion_density *\
                                     np.exp(-self.beta_rad * current_levels['energy'])

    def calculate_tau_sobolev(self, time_exp):
        @np.vectorize
        def get_n_lower(atomic_number, ion_number, level_number):
            level_id = self.levels_dict[(atomic_number, ion_number, level_number)]
            return self.levels[level_id]['number_density']

        n_lower = get_n_lower(self.lines['atomic_number'], self.lines['ion_number'], self.lines['level_id_lower'])

        #(pi*e^2)/(m_e*c)

        sobolev_coefficient = ((np.pi * constants.cgs.e.real ** 2) / (constants.cgs.m_e * constants.cgs.c))
        tau_sobolev = sobolev_coefficient * self.lines['f_lu'] * self.lines['wavelength'] * time_exp * n_lower

        return tau_sobolev


class NebularPlasma(LTEPlasma):
    def __init__(self, abundances, t_rad, w, density, atom_data, density_unit='g/cm^3', t_electron=None):
        Plasma.__init__(self, abundances, t_rad, density, atom_data, density_unit=density_unit)

        self.update_t_rad(t_rad, w, t_electron)

    def update_radiationfield(self, t_rad, w, t_electron=None):
        self.t_rad = t_rad
        self.w = w

        if t_electron is None:
            self.t_electron = 0.9 * self.t_rad

        self.beta_rad = 1 / (t_rad * constants.cgs.k_B)

        self.partition_functions = self._calculate_partition_functions(self.w)

        self.ion_number_density, self.electron_density = self._calculate_ion_populations()


    def calculate_partition_functions(self):
        """
        Calculate partition functions for the ions using the following formula, where
        :math:`i` is the atomic_number, :math:`j` is the ion_number and :math:`k` is the level number.

        .. math::
            Z_{i,j} = \\sum_{k=0}^{max(k)_{i,j}} g_k \\times e^{-E_k / (k_\\textrm{b} T)}



        Returns
        -------

        partition_functions : `~astropy.table.Table`
            with fields atomic_number, ion_number, partition_function

        """

        unique_atom_ion = np.unique(self.levels.__array__()[['atomic_number', 'ion_number']])

        partition_table = []
        for atomic_number, ion_number in unique_atom_ion:
            levels = self.atom_data.get_levels(atomic_number, ion_number)
            partition_function = np.sum(levels['g'][levels['meta']] * np.exp(-levels['energy'][levels['meta']]\
                                                                             * self.beta_rad))
            partition_function += np.sum(levels['g'][~levels['meta']] * np.exp(-levels['energy'][~levels['meta']]\
                                                                               * self.beta_rad))

            partition_table.append((atomic_number, ion_number, partition_function))

        partition_table = np.array(partition_table, dtype=[('atomic_number', np.int),
                                                           ('ion_number', np.int),
                                                           ('partition_function', np.float)])

        partition_table = table.Table(partition_table)
        return partition_table


    def _calculate_phis(self):
        #phis = N(i+1)/N(i) ???
        #calculating ge = 2/(Lambda^3)
        ge = 2 / (np.sqrt(constants.h ** 2 / (2 * np.pi * constants.me * (1 / (self.beta * constants.erg2ev))))) ** 3

        partition_fractions = ge * self.partition_functions[1:] / self.partition_functions[:-1]
        partition_fractions[np.isnan(partition_fractions)] = 0.0
        partition_fractions[np.isinf(partition_fractions)] = 0.0
        #phi = (n_j+1 * ne / nj)
        phis = partition_fractions * np.exp(-self.beta * self.atom_model.ionization_energy)
        zeta = self.atom_model.interpolate_recombination_coefficient(self.t_rad)
        delta = self.atom_model.calculate_radfield_correction_factor(self.t_rad, self.t_electron, self.w)

        phis *= self.w * (delta * zeta + self.w * (1 - zeta)) * (self.t_electron / self.t_rad) ** .5

        return phis

    def calculate_tau_sobolev(self, time_exp):
        wl = self.atom_model.line_list['wl'] * 1e-8
        #C = (np.pi * constants.e**2) / (constants.me * constants.c) #supposed to be (pi*e**2)/(m_e * c)
        Z = self.partition_functions[self.atom_model.line_list['ion'], self.atom_model.line_list['atom'] - 1]
        g_lower = self.atom_model.line_list['g_lower']
        g_upper = self.atom_model.line_list['g_upper']

        e_lower = self.atom_model.line_list['e_lower']
        e_upper = self.atom_model.line_list['e_upper']

        n_lower = np.empty(len(self.atom_model.line_list), dtype=np.float64)
        n_upper = np.empty(len(self.atom_model.line_list), dtype=np.float64)

        ion_number_density = self.ion_number_density[
                             self.atom_model.line_list['ion'], self.atom_model.line_list['atom'] - 1]
        meta = self.atom_model.line_list['metastable']
        non_meta = np.logical_not(meta)
        n_lower[meta] = (g_lower[meta] / Z[meta]) * np.exp(-self.beta * e_lower[meta]) * ion_number_density[meta]
        n_lower[non_meta] = self.w * (g_lower[non_meta] / Z[non_meta]) * np.exp(-self.beta * e_lower[non_meta]) *\
                            ion_number_density[non_meta]

        n_upper[meta] = (g_upper[meta] / Z[meta]) * np.exp(-self.beta * e_upper[meta]) * ion_number_density[meta]
        n_upper[non_meta] = self.w * (g_upper[non_meta] / Z[non_meta]) * np.exp(-self.beta * e_upper[non_meta]) *\
                            ion_number_density[non_meta]
        tau_sobolev_correction = 1 - ((g_lower * n_upper) / (g_upper * n_lower))
        tau_sobolev_correction[np.isnan(tau_sobolev_correction)] = 0.0
        tau_sobolev = constants.sobolev_coeff * self.atom_model.line_list[
                                                'f_lu'] * wl * time_exp * n_lower * tau_sobolev_correction

        return tau_sobolev




