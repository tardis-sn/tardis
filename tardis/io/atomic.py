# atomic model

import os
import logging
import tardis
import numpy as np
import pandas as pd
import cPickle as pickle

from scipy import interpolate
from astropy import units, constants
from collections import OrderedDict
from pandas import DataFrame


class AtomDataNotPreparedError(Exception):
    pass


logger = logging.getLogger(__name__)

tardis_dir = os.path.dirname(os.path.realpath(tardis.__file__))


def data_path(fname):
    return os.path.join(tardis_dir, 'data', fname)


def tests_data_path(fname):
    return os.path.join(tardis_dir, 'tests', 'data', fname)


default_atom_h5_path = data_path('atom_data.h5')


class AtomData(object):
    """
    Class for storing atomic data

    Parameters
    ----------
    atom_data: pandas.DataFrame
        A DataFrame containing the *basic atomic data* with:
            index: atomic_number;
            columns: symbol, name, mass[u].

    ionization_data: pandas.DataFrame
        A DataFrame containing the *ionization data* with:
            index: atomic_number, ion_number;
            columns: ionization_energy[eV].

       It is important to note here is that `ion_number` describes the *final ion state*
            e.g. H I - H II is described with ion=1

    levels_data: pandas.DataFrame
        A DataFrame containing the *levels data* with:
            index: no index;
            columns: atomic_number, ion_number, level_number, energy[eV], g[1], metastable.

    lines_data: pandas.DataFrame
        A DataFrame containing the *lines data* with:
            index: no index;
            columns: line_id, atomic_number, ion_number, level_number_lower, level_number_upper,
                wavelength[angstrom], nu[Hz], f_lu[1], f_ul[1], B_ul[?], B_ul[?], A_ul[1/s].

    macro_atom_data: (pandas.DataFrame, pandas.DataFrame)
        A tuple containing a DataFrame with the *macro atom data* with:
            index: no_index;
            columns: atomic_number, ion_number, source_level_number, destination_level_number,
                transition_line_id, transition_type, transition_probability;
        and a DataFrame with the *macro atom references* with:
            index: no_index;
            columns: atomic_number, ion_number, source_level_number, count_down, count_up, count_total.

        Refer to the docs: http://tardis.readthedocs.io/en/latest/physics/plasma/macroatom.html

    collision_data: (pandas.DataFrame, np.array)
        A tuple containing a DataFrame with the *electron collisions data* with:
            index: atomic_number, ion_number, level_number_lower, level_number_upper;
            columns: e_col_id, delta_e, g_ratio, c_ul;
        and an array with the collision temperatures.

    zeta_data: ?
    synpp_refs: ?
    ion_cx_data: ?

    Attributes
    -------------
    prepared: bool
    has_levels: bool
    has_lines: bool
    has_macro_atom: bool
    atom_data: pandas.DataFrame
    ionization_data: pandas.DataFrame
    macro_atom_data_all: pandas.DataFrame
    macro_atom_references_all: pandas.DataFrame
    collision_data: pandas.DataFrame
    collision_data_temperatures: numpy.array
    symbol2atomic_number: OrderedDict
    atomic_number2symbol OrderedDict

    Methods
    --------
    from_hdfstore
    prepare_atom_data

    Notes
    ------
    1. The units of some columns are given in the square brackets. They are **NOT** the parts of columns' names!

    """

    @classmethod
    def from_hdfstore(cls, fname=None):
        """
        Function to read all the atom data from the special Carsus HDFStore file.

        Parameters
        ----------

        fname: str, optional
            The path to the HDFStore file. If set to `None` the default file with limited atomic_data
            shipped with TARDIS will be used. For more complex atomic data please contact the authors.
            (default: None)
        """

        if fname is None:
            fname = default_atom_h5_path

        if not os.path.exists(fname):
            raise ValueError("Supplied Atomic Model Database %s does not exists" % fname)

        with pd.HDFStore(fname) as store:

            try:
                basic_atom_data = store["basic_atom_df"]
            except KeyError:
                print "Basic atom data is not available in this HDF5-data file."
                basic_atom_data = None

            try:
                ionization_data = store["ionization_df"]
            except KeyError:
                print "Ionization data is not available in this HDF5-data file."
                ionization_data = None

            try:
                levels_data = store["levels_df"]
            except KeyError:
                print "Levels data is not available in this HDF5-data file."
                levels_data = None

            try:
                lines_data = store["lines_df"]
            except KeyError:
                print "Lines data is not available in this HDF5-data file"
                lines_data = None

            try:
                macro_atom_data = store["macro_atom_df"]
            except KeyError:
                print "Macro atom data is not available in this HDF5-data file."
                macro_atom_data = None

            try:
                macro_atom_ref = store["macro_atom_ref_df"]
            except KeyError:
                print "Macro atom reference data is not available in this HDF5-data file."
                macro_atom_ref = None

            try:
                zeta_data = store["zeta_data"]
            except KeyError:
                print "Zeta data is not available in this HDF5-data file."
                zeta_data = None

            try:
                collision_data = store["collisions_df"]
                collision_temperatures = store.get_storer("collisions_df").attrs["temperatures"]
            except KeyError:
                print "Collision data is not available in this HDF5-data file."
                collision_data, collision_temperatures = (None, None)

            try:
                synpp_refs = store["synpp_refs"]
            except KeyError:
                print "Synpp refs is not available in this HDF5-data file."
                synpp_refs = None

            try:
                ion_cx_data = store["ion_cx_data"]
            except KeyError:
                print "Ionization cx data is not available in this HDF5-data file."
                ion_cx_data = None

            atom_data = cls(atom_data=basic_atom_data, ionization_data=ionization_data, levels_data=levels_data,
                            lines_data=lines_data, macro_atom_data=(macro_atom_data, macro_atom_ref), zeta_data=zeta_data,
                            collision_data=(collision_data, collision_temperatures), synpp_refs=synpp_refs,
                            ion_cx_data=ion_cx_data)

            atom_data.uuid1 = store.root._v_attrs['uuid1']
            atom_data.md5 = store.root._v_attrs['md5']

            try:
                atom_data.version = store.root._v_attrs['database_version']
            except KeyError:
                atom_data.version = None

            # ToDo: strore data sources as attributes in carsus
            # if atom_data.version is not None:
            #     atom_data.data_sources = pickle.loads(h5_file.attrs['data_sources'])

            logger.info('Read Atom Data with UUID=%s and MD5=%s', atom_data.uuid1, atom_data.md5)

        return atom_data

    def __init__(self, atom_data, ionization_data, levels_data, lines_data, macro_atom_data=None, zeta_data=None,
                 collision_data=None, synpp_refs=None, ion_cx_data=None):

        self.prepared = False

        self.atom_data = atom_data
        # We have to use constants.u because astropy uses different values for the unit u and the constant.
        # This is changed in later versions of astropy (the value of constants.u is used in all cases)
        if units.u.cgs == constants.u.cgs:
            self.atom_data["mass"] = units.Quantity(self.atom_data["mass"].values, "u").cgs
        else:
            self.atom_data["mass"] = self.atom_data["mass"].values * constants.u.cgs

        self.ionization_data = ionization_data
        self.ionization_data["ionization_energy"] = units.Quantity(self.ionization_data["ionization_energy"].values, "eV").cgs

        if levels_data is not None:
            self.levels = levels_data
            self.levels["energy"] = units.Quantity(self.levels["energy"].values, 'eV').cgs
            self.has_levels = True
        else:
            self.levels = None
            self.has_levels = False

        if lines_data is not None:
            self.lines = lines_data
            self.lines['wavelength_cm'] = units.Quantity(self.lines['wavelength'], 'angstrom').cgs
            self.has_lines = True
        else:
            self.lines = None
            self.has_lines = False

        if macro_atom_data is not None:
            self.macro_atom_data_all, self.macro_atom_references_all = macro_atom_data
            self.has_macro_atom = True
        else:
            self.macro_atom_data_all = None
            self.macro_atom_references_all = None
            self.has_macro_atom = False

        # if ion_cx_data is not None:
        #     self.has_ion_cx_data = True
        #     #TODO:Farm a panda here
        #     self.ion_cx_th_data = DataFrame(np.array(ion_cx_data[0]))
        #     self.ion_cx_th_data.set_index(['atomic_number', 'ion_number', 'level_id'], inplace=True)
        #
        #     self.ion_cx_sp_data = DataFrame(np.array(ion_cx_data[1]))
        #     self.ion_cx_sp_data.set_index(['atomic_number', 'ion_number', 'level_id'])
        # else:
        #     self.has_ion_cx_data = False

        # if zeta_data is not None:
        #     self.zeta_data = zeta_data
        #     self.has_zeta_data = True
        # else:
        #     self.has_zeta_data = False

        if collision_data[0] is not None:
            self.collision_data, self.collision_data_temperatures = collision_data
            self.has_collision_data = True
        else:
            self.collision_data = None
            self.collision_data_temperatures = None
            self.has_collision_data = False

        # if synpp_refs is not None:
        #     self.has_synpp_refs = True
        #     self.synpp_refs = pd.DataFrame(synpp_refs)
        #     self.synpp_refs.set_index(['atomic_number', 'ion_number'], inplace=True)
        #
        # else:
        #     self.has_synpp_refs = False

        self.symbol2atomic_number = OrderedDict(zip(self.atom_data['symbol'].values, self.atom_data.index))
        self.atomic_number2symbol = OrderedDict(zip(self.atom_data.index, self.atom_data['symbol']))

    def prepare_atom_data(self, selected_atomic_numbers, line_interaction_type='scatter', max_ion_number=None,
                          nlte_species=[]):
        """
        Prepares the atom data to set the lines, levels and if requested macro atom data.
        This function mainly cuts the `levels` and `lines` by discarding any data that is not needed (any data
        for atoms that are not needed

        Parameters
        ----------

        selected_atoms : `~set`
            set of selected atom numbers, e.g. set([14, 26])

        line_interaction_type : `~str`
            can be 'scatter', 'downbranch' or 'macroatom'

        max_ion_number : `~int`
            maximum ion number to be included in the calculation

        """
        if not self.prepared:
            self.prepared = True
        else:
            raise AtomDataNotPreparedError("AtomData was already prepared")
        self.selected_atomic_numbers = selected_atomic_numbers

        self.nlte_species = nlte_species
        self.levels = self.levels.reset_index(drop=True)

        self.levels = self.levels[self.levels['atomic_number'].isin(self.selected_atomic_numbers)]

        if max_ion_number is not None:
            self.levels = self.levels[self.levels['ion_number'] <= max_ion_number]

        self.levels = self.levels.set_index(['atomic_number', 'ion_number', 'level_number'])


        self.levels_index = pd.Series(np.arange(len(self.levels), dtype=int), index=self.levels.index)
        #cutting levels_lines
        self.lines = self.lines[self.lines['atomic_number'].isin(self.selected_atomic_numbers)]
        if max_ion_number is not None:
            self.lines = self.lines[self.lines['ion_number'] <= max_ion_number]

        # self.lines.sort(['wavelength', 'line_id'], inplace=True)
        self.lines.sort(['wavelength'], inplace=True)
        self.lines.set_index('line_id', inplace=True)



        self.lines_index = pd.Series(np.arange(len(self.lines), dtype=int), index=self.lines.index)

        tmp_lines_lower2level_idx = pd.MultiIndex.from_arrays([self.lines['atomic_number'], self.lines['ion_number'],
                                                               self.lines['level_number_lower']])

        self.lines_lower2level_idx = self.levels_index.ix[tmp_lines_lower2level_idx].values.astype(np.int64)

        tmp_lines_upper2level_idx = pd.MultiIndex.from_arrays([self.lines['atomic_number'], self.lines['ion_number'],
                                                               self.lines['level_number_upper']])

        self.lines_upper2level_idx = self.levels_index.ix[tmp_lines_upper2level_idx].values.astype(np.int64)

        self.atom_ion_index = None
        self.levels_index2atom_ion_index = None

        if self.has_macro_atom and not (line_interaction_type == 'scatter'):
            self.macro_atom_data = self.macro_atom_data_all[
                self.macro_atom_data_all['atomic_number'].isin(self.selected_atomic_numbers)]

            if max_ion_number is not None:
                self.macro_atom_data = self.macro_atom_data[self.macro_atom_data['ion_number'] <= max_ion_number]

            self.macro_atom_references = self.macro_atom_references_all[
                self.macro_atom_references_all['atomic_number'].isin(
                    self.selected_atomic_numbers)]
            if max_ion_number is not None:
                self.macro_atom_references = self.macro_atom_references[
                    self.macro_atom_references['ion_number'] <= max_ion_number]

            if line_interaction_type == 'downbranch':
                self.macro_atom_data = self.macro_atom_data[(self.macro_atom_data['transition_type'] == -1).values]

                self.macro_atom_references = self.macro_atom_references[self.macro_atom_references['count_down'] > 0]
                self.macro_atom_references['count_total'] = self.macro_atom_references['count_down']
                self.macro_atom_references['block_references'] = np.hstack((0,
                                                                            np.cumsum(self.macro_atom_references[
                                                                                          'count_down'].values[:-1])))
            elif line_interaction_type == 'macroatom':
                block_references = np.hstack((0, np.cumsum(
                    self.macro_atom_references['count_total'].values[:-1])))
                self.macro_atom_references.insert(len(
                    self.macro_atom_references.columns), 'block_references',
                    pd.Series(block_references,
                    index=self.macro_atom_references.index))

            self.macro_atom_references.set_index(['atomic_number', 'ion_number', 'source_level_number'], inplace=True)
            self.macro_atom_references.insert(len(
                    self.macro_atom_references.columns), 'references_idx',
                    pd.Series(np.arange(len(self.macro_atom_references)),
                    index=self.macro_atom_references.index))

            self.macro_atom_data.insert(len(
                self.macro_atom_data.columns), 'lines_idx',
                pd.Series(self.lines_index.ix[self.macro_atom_data[
                'transition_line_id']].values,
                index=self.macro_atom_data.index))

            tmp_lines_upper2level_idx = pd.MultiIndex.from_arrays(
                [self.lines['atomic_number'], self.lines['ion_number'],
                 self.lines['level_number_upper']])

            self.lines_upper2macro_reference_idx = self.macro_atom_references['references_idx'].ix[
                tmp_lines_upper2level_idx].values.astype(np.int64)

            tmp_macro_destination_level_idx = pd.MultiIndex.from_arrays([self.macro_atom_data['atomic_number'],
                                                                         self.macro_atom_data['ion_number'],
                                                                         self.macro_atom_data[
                                                                             'destination_level_number']])

            if line_interaction_type == 'macroatom':
                #Sets all

                self.macro_atom_data.insert(len(
                    self.macro_atom_data.columns), 'destination_level_idx',
                    pd.Series(self.macro_atom_references['references_idx'].ix[
                    tmp_macro_destination_level_idx].values.astype(
                        np.int64), index=self.macro_atom_data.index))

            elif line_interaction_type == 'downbranch':
                # Sets all the destination levels to -1 to indicate that they
                # are not used in downbranch calculations
                self.macro_atom_data.loc[:, 'destination_level_idx'] = (
                    np.ones(len(self.macro_atom_data)) * -1).astype(np.int64)

        self.nlte_data = NLTEData(self, nlte_species)

    def __repr__(self):
        return "<Atomic Data UUID=%s MD5=%s Lines=%d Levels=%d>" % \
               (self.uuid1, self.md5, self.lines.atomic_number.count(), self.levels.energy.count())


class NLTEData(object):
    def __init__(self, atom_data, nlte_species):
        self.atom_data = atom_data
        self.lines = atom_data.lines.reset_index(drop=True)
        self.nlte_species = nlte_species

        if nlte_species:
            logger.info('Preparing the NLTE data')
            self._init_indices()
            self._create_nlte_mask()
            if atom_data.has_collision_data:
                self._create_collision_coefficient_matrix()
        else:
            self._create_nlte_mask()

    def _init_indices(self):
        self.lines_idx = {}
        self.lines_level_number_lower = {}
        self.lines_level_number_upper = {}
        self.A_uls = {}
        self.B_uls = {}
        self.B_lus = {}

        for species in self.nlte_species:
            lines_idx = np.where((self.lines.atomic_number == species[0]) &
                                 (self.lines.ion_number == species[1]))
            self.lines_idx[species] = lines_idx
            self.lines_level_number_lower[species] = self.lines.level_number_lower.values[lines_idx].astype(int)
            self.lines_level_number_upper[species] = self.lines.level_number_upper.values[lines_idx].astype(int)

            self.A_uls[species] = self.atom_data.lines.A_ul.values[lines_idx]
            self.B_uls[species] = self.atom_data.lines.B_ul.values[lines_idx]
            self.B_lus[species] = self.atom_data.lines.B_lu.values[lines_idx]

    def _create_nlte_mask(self):
        self.nlte_levels_mask = np.zeros(self.atom_data.levels.energy.count()).astype(bool)
        self.nlte_lines_mask = np.zeros(self.atom_data.lines.wavelength.count()).astype(bool)

        for species in self.nlte_species:
            current_levels_mask = (self.atom_data.levels.index.get_level_values(0) == species[0]) & \
                           (self.atom_data.levels.index.get_level_values(1) == species[1])
            current_lines_mask = (self.atom_data.lines.atomic_number.values == species[0]) & \
                           (self.atom_data.lines.ion_number.values == species[1])
            self.nlte_levels_mask |= current_levels_mask
            self.nlte_lines_mask |= current_lines_mask


    def _create_collision_coefficient_matrix(self):
        self.C_ul_interpolator = {}
        self.delta_E_matrices = {}
        self.g_ratio_matrices = {}
        collision_group = self.atom_data.collision_data.groupby(level=['atomic_number', 'ion_number'])
        for species in self.nlte_species:
            no_of_levels = self.atom_data.levels.ix[species].energy.count()
            C_ul_matrix = np.zeros((no_of_levels, no_of_levels, len(self.atom_data.collision_data_temperatures)))
            delta_E_matrix = np.zeros((no_of_levels, no_of_levels))
            g_ratio_matrix = np.zeros((no_of_levels, no_of_levels))

            for (atomic_number, ion_number, level_number_lower, level_number_upper), line in \
                collision_group.get_group(species).iterrows():
                C_ul_matrix[level_number_lower, level_number_upper, :] = line.values[2:]
                delta_E_matrix[level_number_lower, level_number_upper] = line['delta_e']
                #TODO TARDISATOMIC fix change the g_ratio to be the otherway round - I flip them now here.
                g_ratio_matrix[level_number_lower, level_number_upper] = line['g_ratio']
            self.C_ul_interpolator[species] = interpolate.interp1d(self.atom_data.collision_data_temperatures,
                                                                   C_ul_matrix)
            self.delta_E_matrices[species] = delta_E_matrix

            self.g_ratio_matrices[species] = g_ratio_matrix


    def get_collision_matrix(self, species, t_electrons):
        c_ul_matrix = self.C_ul_interpolator[species](t_electrons)
        no_of_levels = c_ul_matrix.shape[0]
        c_ul_matrix[np.isnan(c_ul_matrix)] = 0.0

        #TODO in tardisatomic the g_ratio is the other way round - here I'll flip it in prepare_collision matrix

        c_lu_matrix = c_ul_matrix * np.exp(-self.delta_E_matrices[species].reshape((no_of_levels, no_of_levels, 1)) /
                                           t_electrons.reshape((1, 1, t_electrons.shape[0]))) * \
                      self.g_ratio_matrices[species].reshape((no_of_levels, no_of_levels, 1))
        return c_ul_matrix + c_lu_matrix.transpose(1, 0, 2)

