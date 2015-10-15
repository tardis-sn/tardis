import numpy as np
import logging
import os
from pandas import DataFrame
from astropy import units
from astropy import constants as const
import pandas as pd
from tardis.bound_free.exceptions import IncompletePhotoionizationDataError
from scipy.integrate import simps
from tardis import macro_atom

default_photoionization_h5_path = os.path.join(os.path.dirname(__file__), 'data', 'my_phot_data.h5')

logger = logging.getLogger(__name__)


class BaseContinuumData(object):
    def __init__(self, atom_data, photo_dat_fname=None):
        # TODO: remove unnecessary attributes
        self.atom_data = atom_data
        self.levels = atom_data.levels.reset_index()
        self.levels = self.levels.query('atomic_number != ion_number')
        self.continuum_references = self._create_continuum_references()
        self.continuum_data = self._create_continuum_data_from_levels()
        self.macro_atom_data = self._create_macro_atom_data()
        self.photoionization_data = PhotoionizationData.from_hdf5(fname=photo_dat_fname, atom_data=self)
        self.no_levels_with_photdata = len(np.unique(self.photoionization_data.index.values))
        self.multi_index_nu_sorted = self.continuum_data.sort('nu', ascending=False).set_index(
            ['atomic_number', 'ion_number', 'level_number_lower']).index.values
        self.level_number_density = None

    def _create_continuum_data_from_levels(self):
        logger.info('Generating Continuum data from levels data.')
        continuum_data = self.levels.copy(deep=True)
        continuum_data.rename(columns={'level_number': 'level_number_lower', 'index': 'level_lower_index'},
                              inplace=True)
        ionization_data_index = pd.MultiIndex.from_arrays([continuum_data['atomic_number'],
                                                           continuum_data['ion_number'] + 1])
        continuum_data['nu'] = (self.atom_data.ionization_data.ix[ionization_data_index].values.flatten()
                                - continuum_data['energy'].values) * \
                               units.Unit('erg').to('Hz', equivalencies=units.spectral())

        tmp_level_lower_idx = pd.MultiIndex.from_arrays(
            [continuum_data['atomic_number'], continuum_data['ion_number'],
             continuum_data['level_number_lower']])

        continuum_data.insert(len(continuum_data.columns), 'level_lower_idx',
                              pd.Series(self.atom_data.macro_atom_references['references_idx'].ix[
                                  tmp_level_lower_idx].values.astype(
                                  np.int64), index=continuum_data.index))

        tmp_references_idx_index = pd.MultiIndex.from_arrays([continuum_data['atomic_number'],
                                                              continuum_data['ion_number']])

        continuum_data.insert(len(continuum_data.columns), 'continuum_references_idx',
                              pd.Series(self.continuum_references['references_idx'].ix[
                                  tmp_references_idx_index].values.astype(
                                  np.int64), index=continuum_data.index))

        # TODO: this is problematic if there are very close continuum edges
        continuum_data.insert(len(continuum_data.columns), 'continuum_edge_idx',
                              np.arange(len(continuum_data.nu))
                              [(continuum_data.nu.argsort()[::-1]).argsort().values])

        continuum_data.drop(['energy', 'g', 'metastable'], axis=1, inplace=True)

        # continuum_data.reset_index(inplace = True)
        # continuum_data.sort('nu', ascending = False, inplace = True)
        return continuum_data

    def _create_continuum_references(self):
        continuum_references = pd.DataFrame({'counts_total':
                                                 self.levels.reset_index(drop=True).groupby(
                                                     ['atomic_number', 'ion_number']).count().ix[:, 0]})
        continuum_references['counts_total'] *= 2
        block_references = np.hstack((0, np.cumsum(continuum_references['counts_total'].values[:-1])))
        continuum_references.insert(len(continuum_references.columns), 'block_references',
                                    pd.Series(block_references, index=continuum_references.index))

        continuum_references.insert(len(continuum_references.columns), 'references_idx',
                                    pd.Series(np.arange(len(continuum_references)),
                                              index=continuum_references.index))
        return continuum_references

    def _create_macro_atom_data(self):
        macro_atom_continuum_data = self.continuum_data.copy()
        macro_atom_continuum_data.insert(3, 'transition_type',
                                         pd.Series(np.zeros(len(macro_atom_continuum_data), dtype=int),
                                                   index=macro_atom_continuum_data.index))
        target_level_index = pd.MultiIndex.from_arrays([macro_atom_continuum_data['atomic_number'],
                                                        macro_atom_continuum_data['ion_number'],
                                                        macro_atom_continuum_data['level_lower_index']])
        target_level_energy = (self.levels.set_index(['atomic_number', 'ion_number', 'level_number']
        ).loc[target_level_index, 'energy'] *
                               units.Unit('erg').to('eV', equivalencies=units.spectral())).values

        macro_atom_continuum_data.insert(4, 'transition_probability', pd.Series(target_level_energy,
                                                                                index=macro_atom_continuum_data.index))
        tmp = macro_atom_continuum_data.copy()
        tmp['transition_type'] = -3
        tmp['transition_probability'] = tmp.nu * units.Unit('Hz').to('eV', equivalencies=units.spectral())
        macro_atom_continuum_data = pd.concat([macro_atom_continuum_data, tmp])
        macro_atom_continuum_data.sort(['atomic_number', 'ion_number'], ascending=[1, 1])
        return macro_atom_continuum_data


class PhotoionizationData(object):
    @classmethod
    def from_hdf5(cls, fname, atom_data):
        if fname is None:
            fname = default_photoionization_h5_path
        photoionization_data = cls.read_photoionization_data(fname)
        cls.has_needed_photoionization_data(atom_data=atom_data, photoionization_data_all=photoionization_data)
        return photoionization_data

    @staticmethod
    def read_photoionization_data(fname):
        try:
            with pd.HDFStore(fname, 'r') as phot_data:
                photoionization_x_sections = phot_data['photoionization_data']
                return photoionization_x_sections
        except IOError, err:
            print(err.errno)
            print(err)
            raise IOError('Cannot import. Error opening the file to read photoionization cross-sections')

    @staticmethod
    def has_needed_photoionization_data(atom_data, photoionization_data_all):
        # TODO: Instead of testing for the existence of photoionization data for all selected levels, generate
        # continuum data only for levels with photoionization data
        level_indices_all = np.unique(atom_data.continuum_data.set_index(
            ['atomic_number', 'ion_number', 'level_number_lower']).index.values)
        levels_with_photdata_indices = np.unique(photoionization_data_all.index.values)
        mask = np.in1d(level_indices_all, levels_with_photdata_indices)
        if not np.all(mask):
            raise IncompletePhotoionizationDataError(needed_data=level_indices_all[np.logical_not(mask)],
                                                     provided_data=levels_with_photdata_indices,
                                                     list_data_mismatch=True)


class ReturnMCDataMixin(object):
    @property
    def continuum_edges_list(self):
        """
        :return: Array of continuum edge frequencies sorted in decreasing order. Needed for cmontecarlo.
        """
        return self.continuum_data.sort('nu', ascending=False).nu.values

    @property
    def cont_edge2macro_continuum(self):
        """
        :return: Array of continuum_references_idx. Needed for cmontecarlo, to know which macro_atom_continuum
            level corresponds to a certain continuum_edge.
        """
        return self.continuum_data.sort('nu', ascending=False).continuum_references_idx.values

    def get_phot_table_xsect(self, index_nu_sorted):
        # multi_index = self.continuum_data.sort('nu', ascending=False).set_index(
        #    ['atomic_number', 'ion_number', 'level_number_lower']).index.values[index_nu_sorted]
        multi_index = self.multi_index_nu_sorted[index_nu_sorted]
        return self.photoionization_data.loc[multi_index, 'x_sect'].values

    def get_phot_table_nu(self, index_nu_sorted):
        # multi_index = self.continuum_data.sort('nu', ascending=False).set_index(
        #    ['atomic_number', 'ion_number', 'level_number_lower']).index.values[index_nu_sorted]
        multi_index = self.multi_index_nu_sorted[index_nu_sorted]
        return self.photoionization_data.loc[multi_index, 'nu'].values

    def set_level_number_density(self, level_number_density):
        # ? Copy needed
        self.level_number_density = level_number_density.copy().loc[self.multi_index_nu_sorted].values.transpose()


class ContinuumData(BaseContinuumData, ReturnMCDataMixin):
    pass


class TransitionProbabilitiesContinuum(object):
    """
    Outputs:
        transition_probabilities : Pandas DataFrame
    """

    def __init__(self, t_rads, macro_atom_continuum_data,
                 photoionization_data, continuum_references, lte_level_population):
        self.t_rads = t_rads
        self.no_of_shells = len(t_rads)
        self.macro_atom_continuum_data = macro_atom_continuum_data
        self.block_references = np.hstack(
            (continuum_references.block_references, len(self.macro_atom_continuum_data)))
        self.photoionization_data = photoionization_data.copy()
        self.lte_level_population = lte_level_population
        self.data = self.calculate()

    # TODO: maybe resort transition probabilities by value so that less summations are needed in the macro atom
    def calculate(self):
        transition_probabilities = np.zeros((len(self.t_rads), len(self.macro_atom_continuum_data)))
        for i, t_rad in enumerate(self.t_rads):
            transition_probabilities[i] = self._calculate_one_shell(t_rad, self.lte_level_population.loc[:, i])
        transition_probabilities = transition_probabilities.T
        macro_atom.normalize_transition_probabilities(transition_probabilities, self.block_references)
        transition_probabilities = pd.DataFrame(transition_probabilities)
        transition_probabilities.insert(0, 'destination_level_idx',
                                        self.macro_atom_continuum_data.level_lower_idx.values)
        transition_probabilities.insert(1, 'continuum_edge_idx',
                                        self.macro_atom_continuum_data.continuum_edge_idx.values)
        transition_probabilities.insert(2, 'transition_type',
                                        self.macro_atom_continuum_data.transition_type.values)
        return transition_probabilities

    @staticmethod
    def _calculate_one_probability(nu, x_sect, t_rad):
        integrand = (nu ** 2) * x_sect * np.exp(-const.h.cgs.value * nu / (const.k_B.cgs.value * t_rad))
        return simps(integrand, nu)

    def _calculate_one_shell(self, t_rad, lte_levelpop_one_shell):
        transition_probability_row = self.macro_atom_continuum_data['transition_probability'].values.copy()
        tmp = np.zeros(len(self.macro_atom_continuum_data['transition_probability'].values))
        # TODO: ATM each rate is calculated twice; reset_index earlier
        for i, row in self.macro_atom_continuum_data.reset_index().iterrows():
            level_ijk = tuple(map(int, (row['atomic_number'], row['ion_number'], row['level_number_lower'])))
            nu = self.photoionization_data.loc[level_ijk, 'nu'].values
            x_sect = self.photoionization_data.loc[level_ijk, 'x_sect'].values
            tmp[i] = self._calculate_one_probability(nu, x_sect, t_rad) * lte_levelpop_one_shell.loc[level_ijk]
        # TODO: change this in the preparation of the continuum data
        transition_probability_row *= tmp * units.eV.to(units.erg) * 1e-10
        return transition_probability_row
