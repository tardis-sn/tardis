import logging

import numpy as np
import pandas as pd

from scipy import sparse as sp

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.properties.continuum_processes import \
    get_ground_state_multi_index

__all__ = ['MarkovChainTransProbs', 'MarkovChainIndex']

logger = logging.getLogger(__name__)


class SpMatrixSeriesConverterMixin(object):
    @staticmethod
    def series2matrix(series, idx2reduced_idx):
        """
        Converts a Pandas Series to a sparse matrix and re-indexes it.

        Parameters
        ----------
        series : Pandas Series
            Rates or transition probabilities. Indexed by
            source_level_idx, destination_level_idx.
        idx2reduced_idx: Pandas Series
            Values of (compact) matrix index. Indexed by references_idx.
            Maps the references_idx of a level to the index
            used in the sparse matrix.

        Returns
        -------
        matrix : scipy.sparse.coo.coo_matrix
            Sparse matrix of rates or transition probabilites.
        """
        q_indices = (series.index.get_level_values(0),
                     series.index.get_level_values(1))
        q_indices = (idx2reduced_idx.loc[q_indices[0]].values,
                     idx2reduced_idx.loc[q_indices[1]].values)
        max_idx = idx2reduced_idx.max() + 1
        matrix = sp.coo_matrix((series, q_indices), shape=(max_idx, max_idx))
        return matrix

    @staticmethod
    def matrix2series(matrix, idx2reduced_idx, names=None):
        """
        Converts a sparse matrix to a Pandas Series and indexes it.

        Parameters
        ----------
        matrix : scipy.sparse.coo.coo_matrix
            Sparse matrix of rates or transition probabilites.
        idx2reduced_idx: Pandas Series
            Values of (compact) matrix index. Indexed by references_idx.
            Maps the references_idx of a level to the index
            used in the sparse matrix.

        Returns
        -------
        series : Pandas Series
            Rates or transition probabilities. Indexed by
            source_level_idx, destination_level_idx.
        """
        reduced_idx2idx = pd.Series(idx2reduced_idx.index,
                                    index=idx2reduced_idx)
        matrix = matrix.tocoo()
        index = pd.MultiIndex.from_arrays(
            [reduced_idx2idx.loc[matrix.row],
             reduced_idx2idx.loc[matrix.col]]
        )
        series = pd.Series(matrix.data, index=index)
        if names:
            series.index.names = names
        return series


class MarkovChainIndex(ProcessingPlasmaProperty):
    outputs = ('idx2mkv_idx',)

    def calculate(self, atomic_data, continuum_interaction_species):
        ma_ref = atomic_data.macro_atom_references
        mask = ma_ref.index.droplevel('source_level_number').isin(
            continuum_interaction_species
        )
        mask2 = ma_ref.index.isin(get_ground_state_multi_index(
            continuum_interaction_species)
        )
        mask = np.logical_or(mask, mask2)
        idx = ma_ref[mask].references_idx.values
        idx2mkv_idx = pd.Series(np.arange(len(idx)), index=idx)
        idx2mkv_idx.loc['k'] = idx2mkv_idx.max() + 1
        return idx2mkv_idx


class MarkovChainTransProbs(ProcessingPlasmaProperty,
                            SpMatrixSeriesConverterMixin):
    outputs = ('N', 'R', 'B', 'p_deac')

    def calculate(self, p_rad_bb, p_recomb, p_photo_ion, p_coll,
                  idx2mkv_idx):
        p = pd.concat([p_rad_bb, p_recomb, p_photo_ion, p_coll])
        p = p.groupby(level=[0, 1, 2]).sum()
        p = self.normalize_trans_probs(p)
        p_internal = p.xs(0, level='transition_type')
        p_deac = self.normalize_trans_probs(p.xs(-1, level='transition_type'))

        N = pd.DataFrame(columns=p_internal.columns)
        B = pd.DataFrame(columns=p_internal.columns)
        R = pd.DataFrame(columns=p_internal.columns,
                         index=idx2mkv_idx.index)
        R.index.name = 'source_level_idx'
        for column in p_internal:
            Q = self.series2matrix(p_internal[column], idx2mkv_idx)
            inv_N = sp.identity(Q.shape[0]) - Q
            N1 = sp.linalg.inv(inv_N.tocsc())
            R1 = (1 - np.asarray(Q.sum(axis=1))).flatten()
            B1 = N1.multiply(R1)
            N1 = self.matrix2series(N1, idx2mkv_idx,
                                    names=p_internal.index.names)
            B1 = self.matrix2series(B1, idx2mkv_idx,
                                    names=p_internal.index.names)
            N[column] = N1
            B[column] = B1
            R[column] = R1
        N = N.sort_index()
        B = B.sort_index()
        return N, R, B, p_deac

    @staticmethod
    def normalize_trans_probs(p):
        p_summed = p.groupby(level=0).sum()
        index = p.index.get_level_values('source_level_idx')
        p_norm = p / p_summed.loc[index].values
        p_norm = p_norm.fillna(0.0)
        return p_norm
