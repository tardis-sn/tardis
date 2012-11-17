#Calculations of the Plasma conditions

#import constants
import numpy as np
import logging
from astropy import table, units, constants
from collections import OrderedDict

logger = logging.getLogger(__name__)

Bnu = lambda nu, t: (2 * constants.h * nu ** 3 / constants.c ** 2) * np.exp(
    1 / ((constants.h * nu) / (constants.kb * t)))

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

        #Filtering the levels data
        levels_filter = np.zeros(len(self.atom_data._levels)).astype(np.bool)

        for atomic_number in self.abundances['atomic_number']:
            levels_filter = levels_filter | (self.atom_data._levels['atomic_number'] == atomic_number)

        self.levels = self.atom_data._levels[levels_filter]

        #Filtering the lines data

        lines_filter = np.zeros(len(self.atom_data._lines)).astype(np.bool)

        for atomic_number in self.abundances['atomic_number']:
            lines_filter = lines_filter | (self.atom_data._lines['atomic_number'] == atomic_number)

        self.lines = self.atom_data._lines[lines_filter]

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
        Plasma.__init__(self, abundances, t_rad, density, atom_data, density_unit='g/cm^3')

        self.update_t_rad(t_rad)

    def update_t_rad(self, t_rad):
        """
           This functions updates the radiation temperature `t_rad` and calculates the beta_rad
           Parameters. Then calculating :math:`g_e=2/ \\left(\\sqrt{\\frac{h^2}{2\\pi m_e k_\\textrm{B} T}}\\right)^3`.
           Next will calculate the partition functions, followed by the phis
           (using `calculate_saha`).

           Parameters
           ----------
           t_rad : float


       """
        Plasma.update_t_rad(self, t_rad)

        self.partition_functions = self.calculate_partition_functions()

        self.ge = 2 / (np.sqrt(constants.cgs.h.real ** 2 /\
                               (2 * np.pi * constants.cgs.m_e.real * (1 / (self.beta_rad ))))) ** 3

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

            \\Phi_{i, j} = g_e \\times \\frac{Z_{i, j+1}}{Z_{i, j}} e^{\chi_{j\\rightarrow j+1}/k_\\textrm{B}T}




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

"""
class LTEPlasma(Plasma):
    @classmethod
    def from_model(cls, abundances, density, atom_model):
        return cls(abundances, density, atom_model)


    def __init__(self, abundances, density, atom_model):
        \"""
        ionization_energy
        -------------
        ndarray with row being ion, column atom
        
        \"""
        self.atom_model = atom_model
        self.atom_number_density = self._calculate_atom_number_density(abundances, density)
        self.electron_density = np.sum(self.atom_number_density)


    def update_radiationfield(self, t_rad):
        self.t_rad = t_rad
        self.beta = 1 / (t_rad * constants.kbinev)
        self.partition_functions = self._calculate_partition_functions()
        self.ion_number_density, self.electron_density = self._calculate_ion_populations()


    def calculate_tau_sobolev(self, time_exp):
        wl = self.atom_model.line_list['wl'] * 1e-7
        Z = self.partition_functions[self.atom_model.line_list['ion'], self.atom_model.line_list['atom'] - 1]
        g_lower = self.atom_model.line_list['g_lower']
        e_lower = self.atom_model.line_list['e_lower']
        n_lower = (g_lower / Z) * np.exp(-self.beta * e_lower) *\
                  self.ion_number_density[self.atom_model.line_list['ion'],
                                          self.atom_model.line_list['atom'] - 1]
        tau_sobolev = constants.sobolev_coeff * self.atom_model.line_list['f_lu'] * wl * time_exp * n_lower
        return tau_sobolev

    def _calculate_atom_number_density(self, abundances, density):
        if np.abs(np.sum(abundances) - 1) > 1e-10:
            logger.warn('Abundance fractions don\'t add up to one: %s', np.abs(np.sum(abundances) - 1))
        number_density = (density * abundances) / (self.atom_model.masses * constants.u)

        return number_density

    def _calculate_partition_functions(self):
        z_func = lambda energy, g: np.sum(g * np.exp(-self.beta * energy))
        vector_z_func = np.vectorize(z_func, otypes=[np.float64])
        return vector_z_func(self.atom_model.levels_energy, self.atom_model.levels_g)

    def _calculate_phis(self):
        #calculating ge = 2/(Lambda^3)
        ge = 2 / (np.sqrt(constants.h ** 2 / (2 * np.pi * constants.me * (1 / (self.beta * constants.erg2ev))))) ** 3

        partition_fractions = ge * self.partition_functions[1:] / self.partition_functions[:-1]
        partition_fractions[np.isnan(partition_fractions)] = 0.0
        partition_fractions[np.isinf(partition_fractions)] = 0.0
        #phi = (n_j+1 * ne / nj)
        phis = partition_fractions * np.exp(-self.beta * self.atom_model.ionization_energy)
        return phis

    def _calculate_single_ion_populations(self, phis):
        #N1 is ground state
        #N(fe) = N1 + N2 + .. = N1 + (N2/N1)*N1 + (N3/N2)*(N2/N1)*N1 + ... = N1(1+ N2/N1+...)

        ion_fraction_prod = np.cumprod(phis / self.electron_density, axis=0) # (N2/N1, N3/N2 * N2/N1, ...)
        ion_fraction_sum = 1 + np.sum(ion_fraction_prod, axis=0)
        N1 = self.atom_number_density / ion_fraction_sum
        #Further Ns
        Nn = N1 * ion_fraction_prod
        new_electron_density = np.sum(
            Nn * (np.arange(1, self.atom_model.max_ion).reshape((self.atom_model.max_ion - 1, 1))))
        return np.vstack((N1, Nn)), new_electron_density

    def _calculate_ion_populations(self):
        #partition_functions = calculate_partition_functions(energy_data, g_data, beta)
        phis = self._calculate_phis()

        #first estimate
        electron_density = np.sum(self.atom_number_density)
        old_electron_density = self.electron_density
        while True:
            ion_density, electron_density =\
            self._calculate_single_ion_populations(phis)

            if abs(electron_density / old_electron_density - 1) < 0.05: break

            old_electron_density = 0.5 * (old_electron_density + electron_density)
        return ion_density, electron_density

"""

class NebularPlasma(LTEPlasma):
    def update_radiationfield(self, t_rad, w, t_electron=None):
        self.t_rad = t_rad
        self.w = w

        if t_electron is None:
            self.t_electron = 0.9 * self.t_rad

        self.beta = 1 / (t_rad * constants.kbinev)
        self.partition_functions = self._calculate_partition_functions(self.w)
        self.ion_number_density, self.electron_density = self._calculate_ion_populations()

    def _calculate_partition_functions(self, w):
        def calc_partition(energy, g, meta):
            if np.isscalar(meta): return 0
            non_meta = np.logical_not(meta)
            return np.sum(g[meta] * np.exp(-self.beta * energy[meta])) + w * np.sum(
                g[non_meta] * np.exp(-self.beta * energy[non_meta]))
            #return np.sum(g * np.exp(-self.beta * energy))

        vector_calc_partition = np.vectorize(calc_partition, otypes=[np.float64])
        return vector_calc_partition(self.atom_model.levels_energy, self.atom_model.levels_g,
            self.atom_model.levels_metastable)


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


def read_line_list(conn):
    raw_data = []
    curs = conn.execute(
        'select wl, loggf, g_lower, g_upper, e_lower, e_upper, level_id_lower, level_id_upper, atom, ion from lines')
    for wl, loggf, g_lower, g_upper, e_lower, e_upper, level_id_lower, level_id_upper, atom, ion in curs:
        gf = 10 ** loggf
        f_lu = gf / g_lower
        f_ul = gf / g_upper
        raw_data.append((wl, g_lower, g_upper, f_lu, f_ul, e_lower, e_upper, level_id_lower, level_id_upper, atom, ion))

    line_list = np.array(raw_data, dtype=[('wl', np.float64),
                                          ('g_lower', np.int64), ('g_upper', np.int64),
                                          ('f_lu', np.float64), ('f_ul', np.float64),
                                          ('e_lower', np.float64), ('e_upper', np.float64),
                                          ('level_id_lower', np.int64), ('level_id_upper', np.int64),
                                          ('atom', np.int64), ('ion', np.int64)])

    return line_list


def read_ionize_data_from_db(conn, max_atom=30, max_ion=None):
    if max_ion == None: max_ion = max_atom
    ionize_data = np.zeros((max_ion, max_atom))
    ionize_select_stmt = """select
                    atom, ion, ionize_ev
                from
                    ionization
                where
                    atom <= ?
                and
                    ion <= ?"""

    curs = conn.execute(ionize_select_stmt, (max_atom, max_ion))

    for atom, ion, ionize in curs:
        ionize_data[ion - 1, atom - 1] = ionize

    return ionize_data


def read_level_data(conn, max_atom=30, max_ion=None):
    #Constructing Matrix with atoms columns and ions rows
    #dtype is object and the cells will contain arrays with the energy levels
    if max_ion == None:
        max_ion = max_atom
    level_select_stmt = """select
                atom, ion, energy, g, level_id
            from
                levels
            where
                    atom <= ?
                and
                    ion < ?
            order by
                atom, ion, energy"""

    curs = conn.execute(level_select_stmt, (max_ion, max_atom))
    energy_data = np.zeros((max_ion, max_atom), dtype='object')
    g_data = np.zeros((max_ion, max_atom), dtype='object')

    old_elem = None
    old_ion = None

    for elem, ion, energy, g, levelid in curs:
        if elem == old_elem and ion == old_ion:
            energy_data[ion, elem - 1] = np.append(energy_data[ion, elem - 1], energy)
            g_data[ion, elem - 1] = np.append(g_data[ion, elem - 1], g)

        else:
            old_elem = elem
            old_ion = ion
            energy_data[ion, elem - 1] = np.array([energy])
            g_data[ion, elem - 1] = np.array([g])

    return energy_data, g_data


