import numpy as np
import logging
import os
from pandas import DataFrame
from astropy import units
from astropy import constants as const
import pandas as pd
from tardis.continuum.exceptions import IncompletePhotoionizationDataError

default_photoionization_h5_path = os.path.join(os.path.dirname(__file__), 'data', 'my_phot_data.h5')

logger = logging.getLogger(__name__)


class BaseContinuumData(object):
    def __init__(self, atom_data, photo_dat_fname=None):
        # TODO: remove unnecessary attributes
        self.atom_data = atom_data
        self.levels = self._prepare_levels()
        self.continuum_references = self._create_continuum_references()
        self.continuum_data = self._create_continuum_data_from_levels()
        self.macro_atom_data = self._create_macro_atom_data()
        self.photoionization_data = PhotoionizationData.from_hdf5(fname=photo_dat_fname, atom_data=self)
        self.no_levels_with_photdata = len(np.unique(self.photoionization_data.index.values))
        self.multi_index_nu_sorted = self.continuum_data.sort('nu', ascending=False).set_index(
            ['atomic_number', 'ion_number', 'level_number_lower']).index.values
        self.level_number_density = None
        self._set_montecarlo_data()

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

    def _prepare_levels(self):
        levels = self.atom_data.levels.reset_index()
        return levels.query('atomic_number != ion_number')


class ContinuumProcess(object):
    def __init__(self, input_data):
        self.input = input_data

    def _get_level_energy(self, multi_index):
        return (self.input.levels.loc[multi_index, 'energy']).values

    def _get_lte_level_pop(self, multi_index):
        return (self.input.lte_level_pop.loc[multi_index]).values

    def _get_level_pop(self, multi_index):
        return (self.input.level_pop.loc[multi_index]).values

    def _get_level_idx(self, multi_index):
        return self.input.macro_atom_references.loc[multi_index, 'references_idx'].values

    def _get_continuum_idx(self, multi_index_full):
        ion_number_index = self._get_ion_number_index(multi_index_full)
        return self.input.continuum_references.loc[ion_number_index, 'references_idx'].values

    def _get_ion_number_index(self, multi_index_full):
        atomic_number = multi_index_full.get_level_values(0)
        ion_number = multi_index_full.get_level_values(1)
        return pd.MultiIndex.from_arrays([atomic_number, ion_number])

    @property
    def electron_densities(self):
        return self.input.electron_densities

    @property
    def t_rads(self):
        return self.input.t_rads

    @property
    def ws(self):
        return self.input.ws

    @property
    def t_electrons(self):
        return self.input.t_electrons

    @property
    def ion_number_density(self):
        return self.input.ion_number_density

    @property
    def transition_up_filter(self):
        return (self.input.macro_atom_data.transition_type == 1).values

    @property
    def transition_down_filter(self):
        return (self.input.macro_atom_data.transition_type == 0).values

    @property
    def macro_atom_data(self):
        return self.input.macro_atom_data

    @property
    def c_einstein(self):
        return self.input.c_einstein

    @property
    def c0_reg(self):
        return self.input.c_0_regemorter

    @property
    def nu_i(self):
        return self.input.nu_i

    @property
    def photoionization_data(self):
        return self.input.photoionization_data

    @property
    def no_of_shells(self):
        return len(self.input.t_rads)


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


class MCDataMixin(object):
    def _set_montecarlo_data(self):
        nu_sorted_continuum_data = self.continuum_data.sort('nu', ascending=False)
        self.continuum_edges_list = nu_sorted_continuum_data['nu'].values
        self.cont_edge2macro_continuum = nu_sorted_continuum_data['continuum_references_idx'].values

    def get_phot_table_xsect(self, index_nu_sorted):
        multi_index = self.multi_index_nu_sorted[index_nu_sorted]
        return self.photoionization_data.loc[multi_index, 'x_sect'].values

    def get_phot_table_nu(self, index_nu_sorted):
        multi_index = self.multi_index_nu_sorted[index_nu_sorted]
        return self.photoionization_data.loc[multi_index, 'nu'].values

    def set_level_number_density(self, level_number_density):
        level_number_density_tmp = level_number_density.loc[self.multi_index_nu_sorted].values.transpose()
        self.level_number_density = np.ascontiguousarray(level_number_density_tmp)

    def set_level_number_density_ratio(self, level_number_density, lte_level_number_density):
        level_number_density_tmp = level_number_density.loc[self.multi_index_nu_sorted]
        lte_level_number_density_tmp = lte_level_number_density.loc[self.multi_index_nu_sorted]
        level_number_density_ratio = lte_level_number_density_tmp.divide(level_number_density_tmp)
        level_number_density_ratio = level_number_density_ratio.values.transpose()
        self.level_number_density_ratio = np.ascontiguousarray(level_number_density_ratio)


class ContinuumData(BaseContinuumData, MCDataMixin):
    pass
