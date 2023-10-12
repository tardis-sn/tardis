import logging

import numpy as np
import pandas as pd

from scipy import sparse as sp

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.properties.continuum_processes import (
    get_ground_state_multi_index,
)
from tardis.montecarlo.montecarlo_numba.macro_atom import (
    MacroAtomTransitionType,
)

__all__ = [
    "MarkovChainTransProbs",
    "MarkovChainIndex",
    "MarkovChainTransProbsCollector",
    "NonContinuumTransProbsMask",
    "MonteCarloTransProbs",
]

logger = logging.getLogger(__name__)


def normalize_trans_probs(p):
    """
    Normalize a set of transition probabilities.

    Parameters
    ----------
    p : pandas.DataFrame, dtype float
        Unnormalized transition probabilities. Indexed by
        source_level_idx, destination_level_idx.

    Returns
    -------
    pandas.DataFrame, dtype float
        Normalized transition probabilities: the sum of
        all probabilites with the same source_level_idx sum to one.
        Indexed by source_level_idx, destination_level_idx.
    """
    p_summed = p.groupby(level=0).sum()
    index = p.index.get_level_values("source_level_idx")
    p_norm = p / p_summed.loc[index].values
    p_norm = p_norm.fillna(0.0)
    return p_norm


class SpMatrixSeriesConverterMixin(object):
    @staticmethod
    def series2matrix(series, idx2reduced_idx):
        """
        Convert a Pandas Series to a sparse matrix and re-index it.

        Parameters
        ----------
        series : pandas.Series, dtype float
            Rates or transition probabilities. Indexed by
            source_level_idx, destination_level_idx.
        idx2reduced_idx : pandas.Series
            Values of (compact) matrix index. Indexed by references_idx.
            Maps the references_idx of a level to the index
            used in the sparse matrix.

        Returns
        -------
        scipy.sparse.coo.coo_matrix
            Sparse matrix of rates or transition probabilites.
        """
        q_indices = (
            series.index.get_level_values(0),
            series.index.get_level_values(1),
        )
        q_indices = (
            idx2reduced_idx.loc[q_indices[0]].values,
            idx2reduced_idx.loc[q_indices[1]].values,
        )
        max_idx = idx2reduced_idx.max() + 1
        matrix = sp.coo_matrix((series, q_indices), shape=(max_idx, max_idx))
        return matrix

    @staticmethod
    def matrix2series(matrix, idx2reduced_idx, names=None):
        """
        Convert a sparse matrix to a Pandas Series and index it.

        Parameters
        ----------
        matrix : scipy.sparse.coo.coo_matrix
            Sparse matrix of rates or transition probabilites.
        idx2reduced_idx : pandas.Series
            Values of (compact) matrix index. Indexed by references_idx.
            Maps the references_idx of a level to the index
            used in the sparse matrix.
        names : array-like, optional
            Names of levels in MultiIndex of returned Series.

        Returns
        -------
        pandas.Series
            Rates or transition probabilities. Indexed by
            source_level_idx, destination_level_idx.
        """
        reduced_idx2idx = pd.Series(
            idx2reduced_idx.index, index=idx2reduced_idx
        )
        matrix = matrix.tocoo()
        index = pd.MultiIndex.from_arrays(
            [reduced_idx2idx.loc[matrix.row], reduced_idx2idx.loc[matrix.col]]
        )
        series = pd.Series(matrix.data, index=index)
        if names:
            series.index.names = names
        return series


class MarkovChainIndex(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    idx2mkv_idx : pandas.Series, dtype int
    k_packet_idx : int
        Macro atom level idx corresponding to a k-packet.
    idx2deactivation_idx : pandas.Series, dtype int
    """

    outputs = ("idx2mkv_idx", "k_packet_idx", "idx2deactivation_idx")

    def calculate(self, atomic_data, continuum_interaction_species):
        ma_ref = atomic_data.macro_atom_references
        mask = ma_ref.index.droplevel("source_level_number").isin(
            continuum_interaction_species
        )
        mask2 = ma_ref.index.isin(
            get_ground_state_multi_index(continuum_interaction_species)
        )
        mask = np.logical_or(mask, mask2)
        idx = ma_ref[mask].references_idx.values
        idx2mkv_idx = pd.Series(np.arange(len(idx)), index=idx)
        idx2mkv_idx.loc["k"] = idx2mkv_idx.max() + 1

        k_packet_idx = ma_ref.references_idx.max() + 1

        idx2deactivation_idx = idx2mkv_idx + k_packet_idx + 1
        return idx2mkv_idx, k_packet_idx, idx2deactivation_idx


class NonContinuumTransProbsMask(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    non_continuum_trans_probs_mask : numpy.ndarray, dtype bool
    """

    outputs = ("non_continuum_trans_probs_mask",)

    def calculate(self, atomic_data, continuum_interaction_species):
        # I don't have to remove the ground states of
        # the next higher ionization states of the continuum species
        # since they only contain zero probabilities.
        continuum_trans_probs_mask = atomic_data.macro_atom_data.set_index(
            ["atomic_number", "ion_number"]
        ).index.isin(continuum_interaction_species)
        non_continuum_trans_probs_mask = np.logical_not(
            continuum_trans_probs_mask
        )
        return non_continuum_trans_probs_mask


class MarkovChainTransProbsCollector(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    p_combined : pandas.DataFrame, dtype float
        Combined and normalized transition probabilities.
        Indexed by source_level_idx, destination_level_idx.
    """

    outputs = ("p_combined",)

    def __init__(self, plasma_parent, inputs):
        super().__init__(plasma_parent)
        self.inputs = inputs

    def calculate(self, *args):
        p = pd.concat(args)
        p = p.groupby(level=[0, 1, 2]).sum()
        p = normalize_trans_probs(p)
        return p


class MarkovChainTransProbs(
    ProcessingPlasmaProperty, SpMatrixSeriesConverterMixin
):
    outputs = ("N", "R", "B", "p_deactivation")
    latex_name = ("N", "R", "B", r"p_\textrm{deactivation}")
    """
    Attributes
    ----------
    N : pandas.DataFrame, dtype float
        Fundamental matrix of the Markov-chain macro atom.
        Indexed by source_level_idx, destination_level_idx.
        Expected number of visits to destination_level_idx starting
        from souce_level_idx (before being absorbed).
    R : pandas.DataFrame, dtype float
        Deactivation probabilities of the Markov-chain macro atom.
        Indexed by source_level_idx.
        Probability of deactivation/absorption in source_level_idx.
    B : pandas.DataFrame, dtype float
        Absorbing probabilities of the Markov-chain macro atom.
        Indexed by source_level_idx, destination_level_idx.
        Probability of being absorbed in destination_level_idx when
        starting from source_level_idx.
    p_deactivation : pandas.DataFrame, dtype float
        Redistribution probabilities after deactivation of the Markov-chain
        macro atom. Indexed by source_level_idx, destination_level_idx.
        Probability of an r-packet being emitted in the transition
        (source_level_idx --> destination_level_idx) after deactivation
        in source_level_idx.
    """

    def calculate(self, p_combined, idx2mkv_idx):
        p = p_combined
        p_internal = p.xs(0, level="transition_type")
        p_deactivation = normalize_trans_probs(
            p.xs(-1, level="transition_type")
        )

        N = pd.DataFrame(columns=p_internal.columns)
        B = pd.DataFrame(columns=p_internal.columns)
        R = pd.DataFrame(columns=p_internal.columns, index=idx2mkv_idx.index)
        R.index.name = "source_level_idx"
        for column in p_internal:
            Q = self.series2matrix(p_internal[column], idx2mkv_idx)
            inv_N = sp.identity(Q.shape[0]) - Q
            N1 = sp.linalg.inv(inv_N.tocsc())
            R1 = (1 - np.asarray(Q.sum(axis=1))).flatten()
            B1 = N1.multiply(R1)
            N1 = self.matrix2series(
                N1, idx2mkv_idx, names=p_internal.index.names
            )
            B1 = self.matrix2series(
                B1, idx2mkv_idx, names=p_internal.index.names
            )
            N[column] = N1
            B[column] = B1
            R[column] = R1
        N = N.sort_index()
        B = B.sort_index()
        return N, R, B, p_deactivation


class MonteCarloTransProbs(ProcessingPlasmaProperty):
    outputs = (
        "non_continuum_trans_probs",
        "level_absorption_probs",
        "deactivation_channel_probs",
        "transition_probabilities",
        "macro_block_references",
        "macro_atom_data",
    )
    """
    Attributes
    ----------
    non_continuum_trans_probs
    level_absorption_probs
    deactivation_channel_probs
    transition_probabilities
    macro_block_references
    macro_atom_data
    """

    def calculate(
        self,
        non_markov_transition_probabilities,
        atomic_data,
        non_continuum_trans_probs_mask,
        k_packet_idx,
        idx2deactivation_idx,
        level_idxs2transition_idx,
        p_deactivation,
        cool_rate_fb,
        cool_rate_fb_tot,
        level2continuum_idx,
        B,
    ):
        # Prepare the transition probabilities for the non continuum species
        macro_atom_data = atomic_data.macro_atom_data
        transition_info = macro_atom_data[
            ["lines_idx", "transition_type"]
        ].set_index(non_markov_transition_probabilities.index)
        non_continuum_trans_probs = pd.concat(
            [transition_info, non_markov_transition_probabilities], axis=1
        )
        index = macro_atom_data.set_index(
            ["source_level_idx", "destination_level_idx"]
        ).index
        non_continuum_trans_probs = non_continuum_trans_probs.set_index(index)
        non_continuum_trans_probs = non_continuum_trans_probs[
            non_continuum_trans_probs_mask
        ]

        # Prepare the level absorption probabilities for the continuum species
        level_absorption_probs = B.copy()
        level_absorption_probs.insert(0, "lines_idx", -1)
        level_absorption_probs.insert(0, "transition_type", 3)
        destination_level_idx = level_absorption_probs.index.get_level_values(
            "destination_level_idx"
        )
        source_level_idx = level_absorption_probs.rename(
            index={"k": k_packet_idx}
        ).index.get_level_values("source_level_idx")
        destination_level_idx = idx2deactivation_idx.loc[destination_level_idx]
        absorption_index = pd.MultiIndex.from_arrays(
            [source_level_idx, destination_level_idx], names=B.index.names
        )
        level_absorption_probs.index = absorption_index

        # Prepare the free-bound cooling probabilities
        fb_cooling_probs = (
            cool_rate_fb
            / cool_rate_fb_tot.values
            * p_deactivation.loc[("k"), ("bf")]
        )
        continuum_idx = level2continuum_idx.loc[fb_cooling_probs.index].values
        fb_cooling_probs.index = pd.MultiIndex.from_product(
            [["k"], np.ones(len(fb_cooling_probs), dtype=int) * -1],
            names=p_deactivation.index.names,
        )
        fb_cooling_probs.insert(0, "lines_idx", continuum_idx)
        fb_cooling_probs.insert(
            0,
            "transition_type",
            level_idxs2transition_idx.at[("k", "bf"), "transition_type"],
        )

        # Check if there are two-photon decays
        if "two-photon" in p_deactivation.index.get_level_values(1):
            two_photon_index = p_deactivation[
                p_deactivation.index.get_level_values(1) == "two-photon"
            ].index
            level_idxs2transition_idx_two_photon = pd.DataFrame(
                [[-1, MacroAtomTransitionType.TWO_PHOTON.value]],
                index=two_photon_index,
                columns=level_idxs2transition_idx.columns,
            )
            level_idxs2transition_idx = pd.concat(
                [
                    level_idxs2transition_idx_two_photon,
                    level_idxs2transition_idx,
                ]
            )

        # Prepare the deactivation channel probabilities for the continuum species
        deactivation_channel_probs = p_deactivation.copy()
        deactivation_channel_probs = pd.concat(
            [
                level_idxs2transition_idx.reindex(
                    deactivation_channel_probs.index
                ),
                deactivation_channel_probs,
            ],
            axis=1,
        ).reindex(deactivation_channel_probs.index)

        deactivation_channel_probs = deactivation_channel_probs.drop(
            ("k", "bf")
        )
        deactivation_channel_probs = pd.concat(
            [deactivation_channel_probs, fb_cooling_probs], sort=False
        )

        source_level_idx = deactivation_channel_probs.index.get_level_values(
            "source_level_idx"
        )

        source_level_idx = idx2deactivation_idx.loc[source_level_idx]
        destination_level_idx = np.ones_like(source_level_idx) * -1
        deactivation_index = pd.MultiIndex.from_arrays(
            [source_level_idx, destination_level_idx], names=B.index.names
        )

        deactivation_channel_probs.index = deactivation_index

        # Combine everything
        combined_trans_probs = pd.concat(
            [
                level_absorption_probs,
                deactivation_channel_probs,
                non_continuum_trans_probs,
            ]
        )
        combined_trans_probs = combined_trans_probs.sort_index()

        block_references = (
            combined_trans_probs[0].groupby("source_level_idx").count().cumsum()
        )
        continous_index = np.arange(block_references.index.max() + 1)
        block_references = (
            block_references.reindex(continous_index).ffill().astype(int)
        )  # This is needed because some macro atom levels have no transitions
        block_references = np.pad(block_references, (1, 0), constant_values=0.0)

        macro_atom_info = combined_trans_probs[
            ["transition_type", "lines_idx"]
        ].reset_index()

        combined_trans_probs = combined_trans_probs.drop(
            ["lines_idx", "transition_type"], axis="columns"
        )

        return (
            non_continuum_trans_probs,
            level_absorption_probs,
            deactivation_channel_probs,
            combined_trans_probs,
            block_references,
            macro_atom_info,
        )
