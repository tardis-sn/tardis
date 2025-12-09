import logging
import os

from astropy import units as units
import numpy as np
import pandas as pd

from tardis.iip_plasma.continuum.exceptions import (
    IncompletePhotoionizationDataError,
)
from tardis.iip_plasma.continuum.util import *


default_photoionization_h5_path = os.path.join(
    os.path.dirname(__file__), "data", "photoionization_data_H30_He.h5"
)

logger = logging.getLogger(__name__)


class BaseContinuumData(object):
    def __init__(
        self,
        atom_data,
        photo_dat_fname=None,
        selected_continuum_species=[(1, 0), (2, 0)],
    ):
        self.selected_continuum_species = selected_continuum_species
        self.atom_data = atom_data
        self.levels = self._prepare_levels()
        self.continuum_references = self._create_continuum_references()
        self.continuum_data = self._create_continuum_data_from_levels()
        self.photoionization_data = atom_data.photoionization_data
        self.no_levels_with_photdata = len(
            np.unique(self.photoionization_data.index.values)
        )
        self.multi_index_nu_sorted = self.continuum_data.sort_values(
            "nu", ascending=False
        ).index
        self.level_number_density = None
        self._set_montecarlo_data()

    def _create_continuum_data_from_levels(self):
        logger.info("Generating Continuum data from levels data.")
        continuum_data = self.levels.copy(deep=True)
        # TODO: Rethink use of level_lower_index
        continuum_data.rename(
            columns={
                "level_number": "level_number_lower",
                "index": "level_lower_index",
            },
            inplace=True,
        )
        ionization_data_index = self._get_ion_multi_index(
            continuum_data.index, next_higher=True
        )

        continuum_data["nu"] = (
            self.atom_data.ionization_data.loc[ionization_data_index].values
            - continuum_data["energy"].values
        ) * units.Unit("erg").to("Hz", equivalencies=units.spectral())

        continuum_data.insert(
            len(continuum_data.columns),
            "level_lower_idx",
            self.atom_data.macro_atom_references["references_idx"]
            .loc[continuum_data.index]
            .values.astype(np.int64),
        )

        tmp_references_idx_index = self._get_ion_multi_index(
            continuum_data.index, next_higher=False
        )

        continuum_data.insert(
            len(continuum_data.columns),
            "continuum_references_idx",
            self.continuum_references["references_idx"]
            .loc[tmp_references_idx_index]
            .values.astype(np.int64),
        )

        # WARNING: this is problematic if there are very close continuum edges
        continuum_data.insert(
            len(continuum_data.columns),
            "continuum_edge_idx",
            np.arange(len(continuum_data.nu))[
                (continuum_data.nu.argsort()[::-1]).argsort().values
            ],
        )

        continuum_data.drop(["energy", "g", "metastable"], axis=1, inplace=True)

        return continuum_data

    def _create_continuum_references(self):
        continuum_references = pd.DataFrame(
            {
                "counts_total": self.levels.reset_index()
                .groupby(["atomic_number", "ion_number"])
                .count()
                .iloc[:, 0]
            }
        )
        continuum_references["counts_total"] *= 2
        block_references = np.hstack(
            (0, np.cumsum(continuum_references["counts_total"].values[:-1]))
        )
        continuum_references.insert(
            len(continuum_references.columns),
            "block_references",
            block_references,
        )

        continuum_references.insert(
            len(continuum_references.columns),
            "references_idx",
            np.arange(len(continuum_references)),
        )
        return continuum_references

    def _prepare_levels(self):
        mask = self._get_continuum_species_mask(
            self.atom_data, self.selected_continuum_species
        )
        selected_levels = self.atom_data.levels[mask]
        return selected_levels

    @staticmethod
    def _get_continuum_species_mask(atom_data, selected_continuum_species):
        atomic_number = atom_data.levels.index.get_level_values(0)
        ion_number = atom_data.levels.index.get_level_values(1)
        mask = np.zeros_like(atomic_number, dtype=bool)
        for species in selected_continuum_species:
            species_mask = np.logical_and(
                atomic_number == species[0], ion_number == species[1]
            )
            mask = np.logical_or(mask, species_mask)
        return mask

    @staticmethod
    def _get_ion_multi_index(multi_index_full, next_higher=True):
        atomic_number = multi_index_full.get_level_values(0)
        ion_number = multi_index_full.get_level_values(1)
        if next_higher is True:
            ion_number += 1
        return pd.MultiIndex.from_arrays([atomic_number, ion_number])


class MCDataMixin(object):
    def _set_montecarlo_data(self):
        nu_sorted_continuum_data = self.continuum_data.sort_values(
            "nu", ascending=False
        )
        self.continuum_edges_list = nu_sorted_continuum_data["nu"].values
        self.cont_edge2macro_continuum = nu_sorted_continuum_data[
            "continuum_references_idx"
        ].values
        self.bf_emissivities_nd = len(self.photoionization_data["nu"])
        self.bf_emissivities_block_references_old = np.hstack(
            [
                np.array([0]),
                self.photoionization_data["nu"]
                .groupby(level=[0, 1, 2])
                .count()
                .cumsum()
                .values,
            ]
        )
        block_ref = (
            self.photoionization_data["nu"].groupby(level=[0, 1, 2]).count()
        )
        block_ref = block_ref.loc[self.multi_index_nu_sorted].cumsum().values
        self.bf_emissivities_block_references = np.hstack(
            [np.array([0]), block_ref]
        )

    def get_phot_table_xsect(self, index_nu_sorted):
        multi_index = self.multi_index_nu_sorted.values[index_nu_sorted]
        return self.photoionization_data.loc[multi_index, "x_sect"].values

    def get_phot_table_nu(self, index_nu_sorted):
        multi_index = self.multi_index_nu_sorted.values[index_nu_sorted]
        return self.photoionization_data.loc[multi_index, "nu"].values

    def set_level_number_density(self, level_number_density):
        level_number_density_tmp = level_number_density.loc[
            self.multi_index_nu_sorted
        ].values.transpose()
        self.level_number_density = np.ascontiguousarray(
            level_number_density_tmp
        )

    def set_level_number_density_ratio(self, plasma_array):
        level_number_density = plasma_array.level_number_density
        lte_level_number_density = plasma_array.lte_level_number_density
        ion_number_density = plasma_array.ion_number_density
        lte_ion_number_density = plasma_array.lte_ion_number_density

        level_number_density_tmp = level_number_density.loc[
            self.multi_index_nu_sorted
        ]
        lte_level_number_density_tmp = lte_level_number_density.loc[
            self.multi_index_nu_sorted
        ]
        level_number_density_ratio = lte_level_number_density_tmp.divide(
            level_number_density_tmp
        )

        ion_index = get_ion_multi_index(
            self.multi_index_nu_sorted, next_higher=True
        )
        ion_number_density_tmp = ion_number_density.loc[ion_index]
        lte_ion_number_density_tmp = lte_ion_number_density.loc[ion_index]
        ion_number_density_ratio = ion_number_density_tmp.divide(
            lte_ion_number_density_tmp
        )

        level_number_density_ratio = level_number_density_ratio.multiply(
            ion_number_density_ratio.values
        )
        level_number_density_ratio = (
            level_number_density_ratio.values.transpose()
        )

        self.level_number_density_ratio = np.ascontiguousarray(
            level_number_density_ratio
        )


class ContinuumData(BaseContinuumData, MCDataMixin):
    pass


class PhotoionizationData(object):
    @classmethod
    def from_hdf5(cls, fname, atom_data):
        if fname is None:
            fname = default_photoionization_h5_path
        photoionization_data = cls.read_photoionization_data(fname)
        cls.has_needed_photoionization_data(
            atom_data=atom_data, photoionization_data_all=photoionization_data
        )
        selected_atoms = [
            species[0] for species in atom_data.selected_continuum_species
        ]
        # TODO: Fixme! Have to also use ionization stage here
        photoionization_data = photoionization_data.loc[
            pd.IndexSlice[selected_atoms, :, :], :
        ]
        # photoionization_data.drop([(1,0,i) for i in range(14, 20)], inplace=True)
        photoionization_data.loc[(1, 0, 0), "x_sect"] *= 0.0
        try:
            photoionization_data.loc[(2, 0, 0), "x_sect"] *= 0.0
        except:
            print("Cannot set He ground state phot xsect to 0. No He?")
        return cls.asfortranframe(photoionization_data)

    @staticmethod
    def read_photoionization_data(fname):
        try:
            with pd.HDFStore(fname, "r") as phot_data:
                photoionization_x_sections = phot_data["photoionization_data"]
                return photoionization_x_sections
        except (IOError, err):
            print(err.errno)
            print(err)
            raise IOError(
                "Cannot import. Error opening the file to read photoionization cross-sections"
            )

    @staticmethod
    def has_needed_photoionization_data(atom_data, photoionization_data_all):
        level_indices_all = np.unique(atom_data.continuum_data.index.values)
        level_with_photdata_indices = np.unique(
            photoionization_data_all.index.values
        )
        mask = pd.Series(level_indices_all).isin(level_with_photdata_indices)
        if not np.all(mask):
            raise IncompletePhotoionizationDataError(
                needed_data=level_indices_all[np.logical_not(mask)],
                provided_data=level_with_photdata_indices,
                list_data_mismatch=True,
            )

    @staticmethod
    def asfortranframe(df):
        df = pd.DataFrame(
            np.asfortranarray(df.values), columns=df.columns, index=df.index
        )
        return df
