# atomic model

#TODO revisit import statements and reorganize
from scipy import interpolate
import numpy as np
import logging
import os
import h5py

from astropy import table, units

from collections import OrderedDict

from pandas import DataFrame

import pandas as pd

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


def read_hdf5_data(fname, dset_name):
    """This function reads the dataset (dset_name) from the hdf5 file (fname).
    In addition it uses the attribute 'units' and parses it to the `~astropy.table.Table` constructor.

    Parameters
    ----------

    fname : `str`, optional
        path to atomic.h5 file, if set to None it will read in default data directory

    Returns
    -------

    data : `~astropy.table.Table`
        returns the respective
    """

    if fname is None:
        fname = default_atom_h5_path

    h5_file = h5py.File(fname)
    dataset = h5_file[dset_name]
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

    return data_table


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

    data_table = read_hdf5_data(fname, 'basic_atom_data')
    data_table.columns['mass'].convert_units_to('g')

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

    data_table = read_hdf5_data(fname, 'ionization_data')
    data_table.columns['ionization_energy'].convert_units_to('erg')

    return data_table


def read_levels_data(fname=None):
    """This function reads atomic number, ion number, level_number, energy, g, metastable
    information from hdf5 file.

    Parameters
    ----------

    fname : `str`, optional
        path to atomic.h5 file, if set to None it will read in default data directory

    Returns
    -------

    data : `~astropy.table.Table`
        table with fields z[1], ion[1], level_number, energy, g, metastable
    """

    data_table = read_hdf5_data(fname, 'levels_data')
    data_table.columns['energy'].convert_units_to('erg')
    #data_table.columns['ionization_energy'].convert_units_to('erg')

    return data_table


def read_lines_data(fname=None):
    """
    This function reads the wavelength, atomic number, ion number, f_ul, f_l and level id information
     from hdf5 file

    Parameters
    ----------

    fname : `str`, optional
        path to atomic.h5 file, if set to None it will read in default data directory

    Returns
    -------

    data : `~astropy.table.Table`
        table with fields wavelength, atomic_number, ion_number, f_ul, f_lu, level_id_lower, level_id_upper.
    """

    data_table = read_hdf5_data(fname, 'lines_data')
    #data_table.columns['ionization_energy'].convert_units_to('erg')

    return data_table


def read_zeta_data(fname):
    """
    This function reads the recombination coefficient data from the HDF5 file


    :return:
    """

    if fname is None:
        raise ValueError('fname can not be "None" when trying to use NebularAtom')

    if not os.path.exists(fname):
        raise IOError('HDF5 File doesn\'t exist')

    h5_file = h5py.File(fname)

    if 'zeta_data' not in h5_file.keys():
        raise ValueError('zeta_data not available in this HDF5-data file. It can not be used with NebularAtomData')

    zeta_data = h5_file['zeta_data']
    zeta_interp = {}
    t_rads = zeta_data.attrs['t_rad']
    for line in zeta_data:
        zeta_interp[tuple(map(int, line[:2]))] = interpolate.interp1d(t_rads, line[2:])

    return zeta_interp


def read_macro_atom_data(fname):
    if fname is None:
        raise ValueError('fname can not be "None" when trying to use NebularAtom')

    if not os.path.exists(fname):
        raise IOError('HDF5 File doesn\'t exist')

    h5_file = h5py.File(fname)

    if 'macro_atom_data' not in h5_file.keys():
        raise ValueError('Macro Atom Data (macro_atom_data) is not in this HDF5-data file. '
                         'It is needed for complex line interaction')
    macro_atom_data = h5_file['macro_atom_data']

    macro_atom_counts = h5_file['macro_atom_counts']

    return macro_atom_data, macro_atom_counts


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
        """
        Function to read all the atom data from a special TARDIS HDF5 File.

        Parameters
        ----------

        fname: str, optional
            the default for this is `None` and then it will use the very limited atomic_data shipped with TARDIS
            For more complex atomic data please contact the authors.
        """

        atom_data = read_basic_atom_data(fname)
        ionization_data = read_ionization_data(fname)
        levels_data = read_levels_data(fname)
        lines_data = read_lines_data(fname)

        return cls(atom_data=atom_data, ionization_data=ionization_data, levels_data=levels_data,
            lines_data=lines_data)


    def __init__(self, atom_data, ionization_data, levels_data, lines_data):

        self.atom_data = DataFrame(atom_data.__array__())
        self.atom_data.set_index('atomic_number', inplace=True)

        self.ionization_data = DataFrame(ionization_data.__array__())
        self.ionization_data.set_index(['atomic_number', 'ion_number'], inplace=True)

        self.levels_data = DataFrame(levels_data.__array__())
        tmp_levels_index = pd.MultiIndex.from_arrays(levels_data['atomic_number'], levels_data['ion_number'],
            levels_data['level_number'], names=('atomic_number', 'ion_number', 'level_number'))

        self.levels_index = pd.Series(np.arange(len(levels_data), dtype=int), index=tmp_levels_index)


        self.lines_data = DataFrame(lines_data.__array__())
        self.lines_data['nu'] = units.Unit('angstrom').to('Hz', self.lines_data['wavelength'], units.spectral())
        self.lines_data['wavelength_cm'] = units.Unit('angstrom').to('cm', self.lines_data['wavelength'])
        #tmp_lines_index = pd.MultiIndex.from_arrays(self.lines_data)
        #self.lines_inde

        self.symbol2atomic_number = OrderedDict(zip(self.atom_data['symbol'].values, self.atom_data.index))
        self.atomic_number2symbol = OrderedDict(zip(self.atom_data.index, self.atom_data['symbol']))


class NebularAtomData(AtomData):
    """
    Class for storing atomic data. This is a special class geared towards use with the `tardis.plasma.NebularPlasma`-class.
    It contains a way of storing Recombination coefficients and can calculate Zeta values ??STUART please add more info??


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
        """
        Function to read all the atom data from a special TARDIS HDF5 File.

        Parameters
        ----------

        fname: str, optional
            the default for this is `None` and then it will use the very limited atomic_data shipped with TARDIS
            For more complex atomic data please contact the authors.
        """
        zeta_data = read_zeta_data(fname)
        atom_data = read_basic_atom_data(fname)
        ionization_data = read_ionization_data(fname)
        levels_data = read_levels_data(fname)
        lines_data = read_lines_data(fname)

        macro_atom_data, macro_atom_counts = read_macro_atom_data(fname)

        return cls(atom_data=atom_data, ionization_data=ionization_data, levels_data=levels_data,
            lines_data=lines_data, zeta_data=zeta_data, macro_atom_data=(macro_atom_data, macro_atom_counts))


    def __init__(self, atom_data, ionization_data, levels_data, lines_data, zeta_data, macro_atom_data):
        super(NebularAtomData, self).__init__(atom_data, ionization_data, levels_data, lines_data)
        self.zeta_data = zeta_data
        self.macro_atom_data = DataFrame(macro_atom_data[0].__array__())
        self.macro_atom_counts = DataFrame(macro_atom_data[1].__array__())
