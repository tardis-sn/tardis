# atomic model

#TODO revisit import statements and reorganize
from scipy import interpolate
import line
import constants
import numpy as np
import sqlite3
import logging
import os
import h5py

from astropy import table, units

from collections import OrderedDict


try:
    import sqlparse

    sqlparse_available = True
except ImportError:
    sqlparse_available = False

logger = logging.getLogger(__name__)

default_atom_h5_path = os.path.join(os.path.dirname(__file__), 'data', 'atom_data.h5')


@PendingDeprecationWarning
def read_atomic_data(fname=None):
    return read_basic_atom_data(fname)


def read_basic_atom_data(fname=None):
    """This function reads the atomic number, symbol, and mass from hdf5 file

    Parameters
    ----------

    fname : `str`, optional
        path to atomic.h5 file, if set to None it will read in default data directory

    Returns
    -------

    data : `~astropy.table.Table`
        table with fields z[1], symbol, mass[u]
    """

    if fname is None:
        fname = default_atom_h5_path

    h5_file = h5py.File(fname)
    dataset = h5_file['basic_atom_data']
    data = np.asarray(dataset)
    data_units = dataset.attrs['units']

    data_table = table.Table(data)
    for i, col_unit in enumerate(data_units):
        if col_unit == 'n':
            data_table.columns[i].units = None
        elif col_unit == '1':
            data_table.columns[i].units = units.Unit(1)
        else:
            data_table.columns[i].units = units.Unit(col_unit)

    h5_file.close()

    data_table.columns['mass'].convert_units_to('g')

    #converting mass column from u to g

    return data_table


def read_ionization_data(fname=None):
    """This function reads the atomic number, ion number, and ionization energy from hdf5 file

    Parameters
    ----------

    fname : `str`, optional
        path to atomic.h5 file, if set to None it will read in default data directory

    Returns
    -------

    data : `~astropy.table.Table`
        table with fields z[1], ion[1], ionization_energy[eV]
        .. note:: energy from unionized atoms to once-ionized atoms ion = 1, for once ionized
                  to twice ionized ion=2, etc.
    """

    if fname is None:
        fname = default_atom_h5_path

    h5_file = h5py.File(fname)
    dataset = h5_file['ionization_data']
    data = np.asarray(dataset)
    data_units = dataset.attrs['units']

    data_table = table.Table(data)
    for i, col_unit in enumerate(data_units):
        if col_unit == 'n':
            data_table.columns[i].units = None
        elif col_unit == '1':
            data_table.columns[i].units = units.Unit(1)
        else:
            data_table.columns[i].units = units.Unit(col_unit)

    data_table.columns['ionization_energy'].convert_units_to('erg')
    h5_file.close()

    return data_table


def convert_int_ndarray(sqlite_binary):
    if sqlite_binary == '-1':
        return np.array([], dtype=np.int64)
    else:
        return np.frombuffer(sqlite_binary, dtype=np.int64)


def convert_float_ndarray(sqlite_binary):
    if sqlite_binary == '-1.0':
        return np.array([], dtype=np.float64)
    else:
        return np.frombuffer(sqlite_binary, dtype=np.float64)

sqlite3.register_converter('int_ndarray', convert_int_ndarray)
sqlite3.register_converter('float_ndarray', convert_float_ndarray)


class AtomData(object):
    """
    Class for storing atomic data

    AtomData
    ---------

    Parameters
    ----------

    basic_atom_data: ~astropy.table.Table
        containing the basic atom data: z, symbol, and mass

    ionization_data: ~astropy.table.Table
        containing the ionization data: z, ion, and ionization energy
        ::important to note here is that ion describes the final ion state
            e.g. H I - H II is described with ion=2

    """

    @classmethod
    def from_hdf5(cls, fname=None):
        atom_data = read_basic_atom_data(fname)
        ionization_data = read_ionization_data(fname)

        return cls(basic_atom_data=atom_data, ionization_data=ionization_data, levels=None, lines=None)


    def __init__(self, basic_atom_data, ionization_data, levels, lines):
        self.basic_atom_data = basic_atom_data
        self.ionization_data = ionization_data
        self.levels = levels
        self.lines = lines

        self.symbol2atomic_number = OrderedDict(zip(self.basic_atom_data['symbol'], self.basic_atom_data['z']))
        self.atomic_number2symbol = OrderedDict(zip(self.basic_atom_data['atomic_number'], self.basic_atom_data['symbol']))

        self.atom_mass = OrderedDict(self.basic_atom_data['atomic_number', 'mass'])

        self.ionization_data_index = dict(
            [((line['atomic_number'], line['ion_number']), i) for i, line in enumerate(self.ionization_data)])


    def get_ionization_energy(self, atomic_number, ion):
        return self.ionization_data[self.ionization_data_index[atomic_number, ion]]['ionization_energy']


class KuruczAtomModel(AtomData):
    @classmethod
    def from_db(cls, conn, max_atom=30, max_ion=30):
        logger.info('Reading Kurucz model from database max atom=%d max ion=%d', max_atom, max_ion)
        masses = read_atomic_data()['mass'][:max_atom]

        ionization_data = read_ionization_data()


        #reason for max_ion - 1: in energy level data there's unionized, once-ionized, twice-ionized, ...
        #in ionization_energies, there's only once_ionized, twice_ionized
        ionization_energy = np.zeros((max_ion - 1, max_atom))

        for atom, ion, ion_energy in ionization_data:
            if atom > max_atom or ion >= max_ion:
                continue
            ionization_energy[ion - 1, atom - 1] = ion_energy

        levels_energy, levels_g, levels_metastable = read_kurucz_level_data_fromdb(conn, max_atom, max_ion)

        return cls(masses=masses,
            ionization_energy=ionization_energy,
            levels_energy=levels_energy,
            levels_g=levels_g,
            levels_metastable=levels_metastable,
            max_atom=max_atom,
            max_ion=max_ion)

    def __init__(self,
                 masses=None,
                 ionization_energy=None,
                 levels_energy=None,
                 levels_g=None,
                 levels_metastable=None,
                 max_atom=None,
                 max_ion=None):
        self.masses = masses
        self.ionization_energy = ionization_energy
        self.levels_energy = levels_energy
        self.levels_g = levels_g

        self.levels_metastable = levels_metastable
        if max_atom != 30 or max_ion != 30:
            logger.warn('max_atom and/or max_ion are not 30 (max_atom=%s max_ion=%s)', max_atom, max_ion)
        self.max_atom = max_atom
        self.max_ion = max_ion

    def calculate_radfield_correction_factor(self, t_rad, t_electron, w, departure_coefficient=None,
                                             xi_threshold_species=(1, 19)):
    #factor delta ML 1993
        if departure_coefficient is None:
            departure_coefficient = 1 / float(w)
        delta = np.ones((self.max_ion - 1, self.max_atom))
        xi_threshold = self.ionization_energy[xi_threshold_species]

        #Formula 15 ML 1993
        threshold_filter = (self.ionization_energy <= xi_threshold) & (self.ionization_energy > 0)
        delta[threshold_filter] = (t_electron / (departure_coefficient * w * t_rad)) *\
                                  np.exp((delta[threshold_filter] / (constants.kbinev * t_rad))\
                                         - (delta[threshold_filter] / (constants.kbinev * t_electron)))

        threshold_filter = (self.ionization_energy > xi_threshold) & (self.ionization_energy > 0)
        #Formula 20 ML 1993
        delta[self.ionization_energy > xi_threshold] = 1 -\
                                                       np.exp((delta[threshold_filter] / (constants.kbinev * t_rad))\
                                                              - (xi_threshold / (constants.kbinev * t_rad))) +\
                                                       (t_electron / (departure_coefficient * w * t_rad)) *\
                                                       np.exp((delta[threshold_filter] / (constants.kbinev * t_rad))\
                                                              - (
                                                           delta[threshold_filter] / (constants.kbinev * t_electron)))

        return delta


class KuruczMacroAtomModel(KuruczAtomModel):
    @classmethod
    def from_db(cls, conn, max_atom=30, max_ion=30):
        kurucz_atom_model = KuruczAtomModel.from_db(conn, max_atom=max_atom, max_ion=max_ion)
        kurucz_atom_model.macro_atom = line.SimpleMacroAtomData.fromdb(conn)
        return kurucz_atom_model


class CombinedAtomicModel(KuruczAtomModel):
    """Complex Atomic Model, in addition to reading the Kurucz model.
        In addition it reads the recombination to ground state coefficient zeta
        from the zeta table.
        It also includes reading reading the meta stable levels and transition probabilities

        """

    @classmethod
    def from_db(cls, conn, max_atom=30, max_ion=30):
        logger.info('Reading Kurucz model from database max atom=%d max ion=%d', max_atom, max_ion)

        combined_atomic_model = KuruczAtomModel.from_db(conn, max_atom=max_atom, max_ion=max_ion)
        line_list = line.read_line_list(conn, max_atom=max_atom, max_ion=max_ion)
        macro_atom = line.SimpleMacroAtomData.fromdb(conn)

        #factor zeta ML 1993
        recombination_coefficients_t_rads, recombination_coefficients =\
        read_recombination_coefficients_fromdb(conn,
            max_atom,
            max_ion)

        return cls(masses=combined_atomic_model.masses,
            ionization_energy=combined_atomic_model.ionization_energy,
            levels_energy=combined_atomic_model.levels_energy,
            levels_g=combined_atomic_model.levels_g,
            levels_metastable=combined_atomic_model.levels_metastable,
            line_list=line_list,
            macro_atom=macro_atom,
            recombination_coefficients_t_rads=recombination_coefficients_t_rads,
            recombination_coefficients=recombination_coefficients,
            max_atom=max_atom,
            max_ion=max_ion)


    def __init__(self,
                 masses=None,
                 ionization_energy=None,
                 levels_energy=None,
                 levels_g=None,
                 levels_metastable=None,
                 line_list=None,
                 macro_atom=None,
                 recombination_coefficients_t_rads=None,
                 recombination_coefficients=None,
                 max_atom=None,
                 max_ion=None):
        KuruczAtomModel.__init__(self,
            masses=masses,
            ionization_energy=ionization_energy,
            levels_energy=levels_energy,
            levels_g=levels_g,
            levels_metastable=levels_metastable,
            max_atom=max_atom,
            max_ion=max_ion
        )
        self.recombination_coefficients_t_rads = recombination_coefficients_t_rads
        self.recombination_coefficients = recombination_coefficients
        self.line_list = line_list
        self.macro_atom = macro_atom

    def interpolate_recombination_coefficient(self, t_rad, kind='linear', bounds_error=True, fill_value=1.):
        interpolator = interpolate.interp1d(self.recombination_coefficients_t_rads,
            self.recombination_coefficients,
            kind=kind,
            bounds_error=bounds_error,
            fill_value=fill_value)
        return interpolator(t_rad)


def read_recombination_coefficients_fromdb(conn, max_atom=30, max_ion=30):
    select_zeta_stmt = """select
                            atom, ion, zeta
                        from
                            zeta
                        where
                                atom < %d
                            and
                                ion < %d
                        """ % (max_atom, max_ion)

    logger.debug('Reading recombination coefficients from db:\n%s\n%s\n%s',
        '-' * 80,
        select_zeta_stmt,
        '-' * 80)

    curs = conn.execute(select_zeta_stmt)
    t_rads = np.arange(2000, 42000, 2000)
    recombination_coefficients = np.ones((max_ion - 1, max_atom, len(t_rads)))
    for atom, ion, zeta in curs:
        recombination_coefficients[ion - 1, atom - 1] = zeta
        #interpolator = interpolate.interp1d(t_rads, recombination_coefficients, kind='linear', bounds_error=False,
    #    fill_value=1.)
    return t_rads, recombination_coefficients


def read_kurucz_level_data_fromdb(conn, max_atom=30, max_ion=None):
    #Constructing Matrix with atoms columns and ions rows
    #dtype is object and the cells will contain arrays with the energy levels
    if max_ion == None:
        max_ion = max_atom

    level_select_stmt = """select
                atom, ion, energy, g, metastable, level_id
            from
                levels
            where
                    atom <= %d
                and
                    ion < %d
            order by
                atom, ion, energy""" % (max_atom, max_ion)

    if sqlparse_available:
        logger.debug('Reading level data from db:\n%s\n%s\n%s',
            '-' * 80,
            sqlparse.format(level_select_stmt,
                reindent=True,
                keyword_case='upper'),

            '-' * 80)
    else:
        logger.debug(level_select_stmt)

    curs = conn.execute(level_select_stmt)
    energy_data = np.zeros((max_ion, max_atom), dtype='object')
    g_data = np.zeros((max_ion, max_atom), dtype='object')
    metastable_data = np.zeros((max_ion, max_atom), dtype='object')

    old_elem = None
    old_ion = None

    for elem, ion, energy, g, metastable, levelid in curs:
        if elem == old_elem and ion == old_ion:
            energy_data[ion, elem - 1] = np.append(energy_data[ion, elem - 1], energy)
            g_data[ion, elem - 1] = np.append(g_data[ion, elem - 1], g)
            metastable_data[ion, elem - 1] = np.append(metastable_data[ion, elem - 1], np.bool(metastable))
        else:
            old_elem = elem
            old_ion = ion
            energy_data[ion, elem - 1] = np.array([energy])
            g_data[ion, elem - 1] = np.array([g])
            metastable_data[ion, elem - 1] = np.array([np.bool(metastable)], dtype=np.bool)

    return energy_data, g_data, metastable_data
