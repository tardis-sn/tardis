#Calculations of the Plasma conditions

#import constants
import numpy as np
import logging
from astropy import table, units, constants
from collections import OrderedDict

from pandas import DataFrame, Series, Index
import pandas as pd

logger = logging.getLogger(__name__)
import pdb
#Bnu = lambda nu, t: (2 * constants.h * nu ** 3 / constants.c ** 2) * np.exp(
#    1 / ((constants.h * nu) / (constants.kb * t)))


#Defining soboleve constant

sobolev_coefficient = ((np.pi * constants.cgs.e.value ** 2) / (constants.cgs.m_e.value * constants.cgs.c.value))

class Plasma(object):
    """
    Model for Plasma

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
    #TODO make density a astropy.quantity
    def __init__(self, abundances, t_rad, density, atom_data, density_unit='g/cm^3', max_ion_number=None):
        self.atom_data = atom_data
        self.density = density

        self.abundances = self.calculate_atom_number_densities(abundances)

        # Filtering the levels & lines data for unnecesary atoms
        levels_atom_filter = atom_data.levels_data['atomic_number'].isin(self.abundances.index)
        self.levels_data = atom_data.levels_data[levels_atom_filter]

        lines_atom_filter = atom_data.lines_data['atomic_number'].isin(self.abundances.index)
        self.lines_data = atom_data.lines_data[lines_atom_filter]

        if max_ion_number is not None:
            self.levels_data = self.levels_data[self.levels_data['ion_number'] <= max_ion_number]
            self.lines_data = self.lines_data[self.lines_data['ion_number'] <= max_ion_number]


        # Setting the indices
        self.levels_data.set_index(['atomic_number', 'ion_number', 'level_number'], inplace=True)
        self.lines_data.set_index(['atomic_number', 'ion_number', 'level_number_lower'], inplace=True)


    def validate_atom_data(self):
        required_attributes = ['lines', 'levels']
        for attribute in required_attributes:
            if not hasattr(self.atom_data, attribute):
                raise ValueError('AtomData incomplete missing')


    def calculate_atom_number_densities(self, abundances):
        """
        Calculates the atom number density, using the following formula, where Z is the atomic number
        and X is the abundance fraction

        .. math::
            N_{Z} = \\frac{\\rho_\\textrm{total}\\times \\textrm{X}_\\textrm{Z}}{m_\\textrm{Z}}

        """

        #Converting abundances
        atomic_numbers = [self.atom_data.symbol2atomic_number[key] for key in abundances.keys()]
        abundance_fractions = Series(abundances.values(),
            index=Index(atomic_numbers, dtype=np.int, name='atomic_number'),
            name='abundance_fraction')



        #Normalizing Abundances

        abundance_sum = abundance_fractions.sum()

        if abs(abundance_sum - 1) > 1e-5:
            logger.warn('Abundances do not add up to 1 (Sum = %.4f). Renormalizing', (abundance_sum))

        abundance_fractions /= abundance_sum

        number_densities = (abundance_fractions * self.density) /\
                           self.atom_data.atom_data.ix[abundance_fractions.index]['mass']

        return DataFrame({'abundance_fraction': abundance_fractions, 'number_density': number_densities})

    def update_radiationfield(self, t_rad):
        """
        This functions updates the radiation temperature `t_rad` and calculates the beta_rad
        Parameters
        ----------
        t_rad : float


        """
        self.t_rad = t_rad
        self.beta_rad = 1 / (constants.cgs.k_B.value * t_rad)


class LTEPlasma(Plasma):
    """
    Model for Plasma using a local thermodynamic equilibrium approximation.

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

    def __init__(self, abundances, t_rad, density, atom_data, density_unit='g/cm^3', max_ion_number=None):
        Plasma.__init__(self, abundances, t_rad, density, atom_data, density_unit=density_unit,
            max_ion_number=max_ion_number)
        self.ion_number_density = None
        #self.update_radiationfield(t_rad)


    def update_radiationfield(self, t_rad, n_e_convergence_threshold=0.05):
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
        Plasma.update_radiationfield(self, t_rad)

        self.partition_functions = self.calculate_partition_functions()

        self.ge = ((2 * np.pi * constants.cgs.m_e.value / self.beta_rad) / (constants.cgs.h.value ** 2)) ** 1.5

        #Calculate the Saha ionization balance fractions
        phis = self.calculate_saha()

        #initialize electron density with the sum of number densities
        electron_density = self.abundances['number_density'].sum()
        self.calculate_ionization_balance(phis, electron_density)

        n_e_iterations = 0

        while True:
            self.calculate_ionization_balance(phis, electron_density)
            ion_numbers = np.array([item[1] for item in self.ion_number_density.index])
            new_electron_density = np.sum(self.ion_number_density.values * ion_numbers)
            n_e_iterations += 1
            if abs(new_electron_density - electron_density) / electron_density < n_e_convergence_threshold: break
            electron_density = 0.5 * (new_electron_density + electron_density)

        logger.info('Took %d iterations to converge on electron density' % n_e_iterations)

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


        def group_calculate_partition_function(group):
            return np.sum(group['g'] * np.exp(-group['energy'] * self.beta_rad))

        partition_functions = self.levels_data.groupby(level=['atomic_number', 'ion_number']).apply(
            group_calculate_partition_function)
        return partition_functions

    def calculate_saha(self):
        """
        Calculating the ionization equilibrium using the Saha equation, where i is atomic number,
        j is the ion_number, :math:`n_e` is the electron density, :math:`Z_{i, j}` are the partition functions
        and :math:`\chi` is the ionization energy.

        .. math::


            \\Phi_{i,j} = \\frac{N_{i, j+1} n_e}{N_{i, j}}

            \\Phi_{i, j} = g_e \\times \\frac{Z_{i, j+1}}{Z_{i, j}} e^{-\chi_{j\\rightarrow j+1}/k_\\textrm{B}T}

        """
        phis = {}
        for atomic_number, part_funcs in self.partition_functions.groupby(level=['atomic_number']):
            part_funcs = part_funcs.dropna()

            ionization_energy = self.atom_data.ionization_data.ix[part_funcs.index[1:]]['ionization_energy'].values
            current_phis = self.ge * (part_funcs[1:].values / part_funcs[:-1].values) * np.exp(
                -self.beta_rad * ionization_energy)
            phis[atomic_number] = current_phis
            #pdb.set_trace()

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

        for atomic_number in phis:
            current_phis = phis[atomic_number] / electron_density
            phis_product = np.cumproduct(current_phis)

            neutral_atom_density = self.abundances.ix[atomic_number]['number_density'] / (1 + np.sum(phis_product))
            ion_densities = [neutral_atom_density] + list(neutral_atom_density * phis_product)

            self.ion_number_density.ix[atomic_number] = ion_densities


    def calculate_level_populations(self):
        """
        Calculate the level populations and putting them in the column 'number-density' of the self.levels table.
        :math:`N` denotes the ion number density calculated with `calculate_ionization_balance`, i is the atomic number,
        j is the ion number and k is the level number.

        .. math::
            \\frac{g_k}{Z_{i,j}} \\times N_{i, j} \\times e^{-\\beta_\\textrm{rad} \\times E_k}



        :return:
        """

        atom_ion_index = self.levels_data.index.droplevel(level=2)

        #partition_functions = Z

        Z = self.partition_functions.ix[atom_ion_index]

        ion_number_density = self.ion_number_density.ix[atom_ion_index]

        level_populations = (self.levels_data['g'].values / Z.values) * ion_number_density.values *\
                            np.exp(-self.beta_rad * self.levels_data['energy'].values)

        self.level_populations = Series(level_populations, index=self.levels_data.index)

    def calculate_tau_sobolev(self, time_exp):
        """
        This function calculates the Sobolev optical depth :math:`\\tau_\\textrm{Sobolev}`



        .. math::
            C_\\textrm{Sobolev} = \\frac{\\pi e^2}{m_e c}

            \\tau_\\textrm{Sobolev} = C_\\textrm{Sobolev}\,  \\lambda\\, f_{\\textrm{lower}\\rightarrow\\textrm{upper}}\,
                t_\\textrm{explosion}\, N_\\textrm{lower}



        .. note::
            Currently we're ignoring the term for stimulated emission:
                :math:`(1 - \\frac{g_\\textrm{lower}}{g_\\textrm{upper}}\\frac{N_\\textrm{upper}}{N_\\textrm{lower}})`


        """

        tau_sobolev = sobolev_coefficient * self.lines_data['f_lu'] * self.lines_data['wavelength'] *\
                      time_exp * self.level_populations.ix[self.lines_data.index].values

        return tau_sobolev


class NebularPlasma(LTEPlasma):
    """
    Model for Plasma using the Nebular approximation

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

    def __init__(self, abundances, t_rad, w, density, atom_data, t_electron=None, density_unit='g/cm^3'):
        Plasma.__init__(self, abundances, t_rad, density, atom_data, density_unit=density_unit)
        self.update_radiationfield(t_rad, t_electron, w)


    def update_radiationfield(self, t_rad, t_electron, w, n_e_convergence_threshold=0.05):
        self.t_rad = t_rad
        self.w = w

        if t_electron is None:
            self.t_electron = 0.9 * self.t_rad

        self.beta_rad = 1 / (t_rad * constants.cgs.k_B.value)
        self.beta_electron = 1 / (self.t_electron * constants.cgs.k_B.value)

        self.partition_functions = self.calculate_partition_functions()

        self.ge = ((2 * np.pi * constants.cgs.m_e.value / self.beta_rad) / (constants.cgs.h.value ** 2)) ** 1.5

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


        #self.ion_number_density, self.electron_density = self._calculate_ion_populations()


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

        unique_atom_ion = np.unique(self.levels.__array__()[['atomic_number', 'ion_number']])

        partition_table = []
        for atomic_number, ion_number in unique_atom_ion:
            levels = self.atom_data.get_levels(atomic_number, ion_number)
            partition_function = np.sum(
                levels['g'][levels['metastable']] * np.exp(-levels['energy'][levels['metastable']]\
                                                           * self.beta_rad))
            partition_function += self.w * np.sum(
                levels['g'][~levels['metastable']] * np.exp(-levels['energy'][~levels['metastable']]\
                                                            * self.beta_rad))

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

        phi_table = super(NebularPlasma, self).calculate_saha()
        phi_table = self.calculate_radiation_field_correction(phi_table)
        zetas = []
        for line in phi_table:
            zetas.append(self.atom_data.zeta_data[(line['atomic_number'], line['numerator'])](self.t_rad))

        phi_table.add_column(table.Column('zeta', zetas))

        phi_table['phi'] *= self.w * (phi_table['delta'][:, 0] * phi_table['zeta'] + self.w * (1 - phi_table['zeta'])) *\
                            (self.t_electron / self.t_rad) ** .5

        return phi_table


    def calculate_radiation_field_correction(self, phi_table, departure_coefficient=None,
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

        chi_threshold = self.atom_data.get_ions(*chi_threshold_species)['ionization_energy']
        deltas = []
        for line in phi_table:
            current_chi = self.atom_data.get_ions(line['atomic_number'], line['numerator'])['ionization_energy']

            delta = (self.t_electron / (self.w * self.t_rad)) *\
                    np.exp(self.beta_rad * chi_threshold - self.beta_electron * current_chi)
            #TODO check this Stuart - there seems to be an error in Formula 15 in ML93
            if chi_threshold >= current_chi:
                deltas.append(delta)


            else:
                delta = 1 - np.exp(self.beta_rad * chi_threshold - self.beta_rad * current_chi) + delta
                deltas.append(delta)

        phi_table.add_column(table.Column('delta', deltas))

        return phi_table


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
            current_number_density = n_levels[level_filter]


            #boolean index array
            current_meta = current_levels['metastable']
            current_number_density[current_meta] = (current_levels['g'][current_meta] / partition_function) *\
                                                   current_ion_density *\
                                                   np.exp(-self.beta_rad * current_levels['energy'][current_meta])

            current_number_density[~current_meta] = self.w * (current_levels['g'][~current_meta] / partition_function) *\
                                                    current_ion_density *\
                                                    np.exp(-self.beta_rad * current_levels['energy'][~current_meta])
