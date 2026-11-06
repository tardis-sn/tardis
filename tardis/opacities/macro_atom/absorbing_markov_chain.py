"""Calculate first-exit probabilities for reduced macro-atom chains.

This module eliminates internal transitions for selected blocks independently in each cell
and computes the source-state probabilities for retained transitions.

Effectively this compresses the macro atom evaluation by internal reducing probability chains
to equivalent single appropriate probability steps.

References
----------
    Absorbing Markov chain theory: https://en.wikipedia.org/wiki/Absorbing_Markov_chain
"""

from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import NamedTuple

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import splu

from tardis.util.base import element_symbol2atomic_number


class AbsorbingChainResult(NamedTuple):
    """
    Carry processed reduced block exit probabilities and their global/compact index mapping.

    Necessary for the monte-carlo step to enter and exit the reduced block.
    """

    matrix: npt.NDArray[np.float64]
    residual_probabilities: pd.DataFrame
    macro_level_to_reduced_idx: npt.NDArray[np.int64]
    reduced_to_macro_level_idx: npt.NDArray[np.int64]


def create_absorbing_probs(
    transition_probabilities: pd.DataFrame,
    metadata: pd.DataFrame,
    mode: str = "full",
    elements: Sequence[str] = (),
    include_thermal_pool: bool = False,
    include_ionization_pool: bool = False,
    max_workers: int = 1,
) -> AbsorbingChainResult:
    """Calculate absorbing Markov chain probabilities and deactivation data.

    Computes the absorption probability matrices for each cell and extracts
    deactivation probabilities by solving the absorbing Markov chain system
    for the requested block.
    Eliminates internal transition inside handled blocks.

    Parameters
    ----------
    transition_probabilities : pd.DataFrame
        DataFrame with transition probabilities between states for each cell.
        Rows represent transitions, columns represent cells.
    metadata : pd.DataFrame
        Metadata about transitions including source_level_idx,
        destination_level_idx, and transition_type. Values >= 0 for
        transition_type indicate internal transitions transitions that are not accompanied
        by reemission of the r-packet, and an end to the interaction handler chain.
        Negative values indicate deactivation, or reemitting from an absorbing state.
    mode : {'full', 'partial', 'explicit'}
        Full eliminates all internal jumps. Partial eliminates only jumps with
        both endpoints selected (e.g., internal transitions for a requested element),
        optionally including either shared pool.
        Explicit retains every original transition.
    elements : sequence of str
        Element symbols, required only in partial mode. All represented ions
        and levels of each selected element are included.
    include_thermal_pool : bool, optional
        Include the shared k state in partial reduction. Defaults to False.
    include_ionization_pool : bool, optional
        Include the shared i state in partial reduction. Defaults to False.
    max_workers : int
        Maximum number of cell-independent sparse LU workers.

    Returns
    -------
    AbsorbingChainResult
        Cell-by-activation-by-exit-source matrix, containing deactivation probabilities,
        and explicit global-to-reduced and reduced-to-global int64 maps.

    Notes
    -----
    The solve is ``(I-Q) B = diag(r)``, where r is the sum of actual remaining
    channels. Zero residual rows stay zero.

    The matrix selects the source of the first retained transition, which can
    itself be internal. It is not necessarily a final emission state because
    only requested physical blocks are blocks are processed out of the matrix,
    and the not-requested blocks retain full internal transition probability tables.
    This means you can move freely back and forth between the absorbing-processed block.
    """
    if mode not in ("full", "partial", "explicit"):
        raise ValueError("markov_chain.mode must be full, partial, or explicit")
    if (mode == "partial") != bool(elements):
        raise ValueError(
            "markov_chain.elements must not be empty in partial mode. Specify which elements to process"
        )
    if mode != "partial" and (include_thermal_pool or include_ionization_pool):
        raise ValueError(
            "markov_chain pool selection requires mode='partial'; full mode "
            "already includes both pools and explicit mode reduces neither"
        )
    try:
        atomic_numbers = sorted(
            {element_symbol2atomic_number(symbol) for symbol in elements}
        )
    except Exception as error:
        raise ValueError(
            f"Invalid markov_chain element symbol in {list(elements)}"
        ) from error

    sources = metadata.source_level_idx.to_numpy(dtype=np.int64)
    destinations = metadata.destination_level_idx.to_numpy(dtype=np.int64)
    internal = metadata.transition_type.to_numpy() >= 0
    global_size = int(sources.max()) + 1
    source_states = np.unique(sources)
    missing = np.setdiff1d(destinations[internal], source_states)
    if len(missing):
        raise ValueError(
            f"Internal destination levels lack outgoing transition blocks: {missing}"
        )
    if mode == "full":
        selected = np.arange(global_size, dtype=np.int64)
    elif mode == "partial":
        # Source tuples start with an atomic number or a shared-pool label.
        source_identifiers = metadata.source.str[0]
        selected_rows = source_identifiers.isin(atomic_numbers).to_numpy()
        represented = set(source_identifiers)
        absent = set(atomic_numbers) - represented
        if absent:
            raise ValueError(
                f"markov_chain.elements absent from assembled macro-atom network: {sorted(absent)}"
            )
        if include_thermal_pool:
            selected_rows = (
                selected_rows | (metadata.source == ("k", -99, -99)).to_numpy()
            )
        if include_ionization_pool:
            selected_rows = (
                selected_rows | (metadata.source == ("i", -99, -99)).to_numpy()
            )
        selected = np.unique(sources[selected_rows])
    else:  # Explicit mode skips the calculation and reduction
        selected = np.empty(0, dtype=np.int64)
    lookup = np.full(global_size, -1, dtype=np.int64)
    lookup[selected] = np.arange(len(selected), dtype=np.int64)
    # Negative deactivation destinations must never index the lookup.
    eliminated = np.zeros(len(metadata), dtype=bool)
    eliminated[internal] = (lookup[sources[internal]] >= 0) & (
        lookup[destinations[internal]] >= 0
    )
    residual = transition_probabilities.copy()
    residual.iloc[eliminated, :] = 0.0
    rows = lookup[sources[eliminated]]
    columns = lookup[destinations[eliminated]]
    matrices = np.zeros(
        (transition_probabilities.shape[1], len(selected), len(selected))
    )

    if mode == "explicit":
        return AbsorbingChainResult(matrices, residual, lookup, selected)

    internal_matrices = (
        sparse.coo_matrix(
            (
                transition_probabilities.iloc[:, cell].to_numpy()[eliminated],
                (rows, columns),
            ),
            shape=(len(selected), len(selected)),
        )
        for cell in range(len(matrices))
    )
    if mode == "full":
        if max_workers > 1 and len(matrices) > 1:
            with ThreadPoolExecutor(
                max_workers=min(max_workers, len(matrices))
            ) as executor:
                for cell, solved in enumerate(
                    executor.map(
                        solve_cell_full,
                        internal_matrices,
                        buffersize=max_workers,
                    )
                ):
                    matrices[cell] = solved
        else:
            for cell, solved in enumerate(
                map(solve_cell_full, internal_matrices)
            ):
                matrices[cell] = solved
    elif mode == "partial":
        exit_rates = (
            np.bincount(
                sources,
                weights=residual.iloc[:, cell].to_numpy(),
                minlength=global_size,
            )[selected]
            for cell in range(len(matrices))
        )
        if max_workers > 1 and len(matrices) > 1:
            with ThreadPoolExecutor(
                max_workers=min(max_workers, len(matrices))
            ) as executor:
                for cell, solved in enumerate(
                    executor.map(
                        solve_cell_partial,
                        internal_matrices,
                        exit_rates,
                        buffersize=max_workers,
                    )
                ):
                    matrices[cell] = solved
        else:
            for cell, solved in enumerate(
                map(
                    solve_cell_partial,
                    internal_matrices,
                    exit_rates,
                    strict=True,
                )
            ):
                matrices[cell] = solved
    return AbsorbingChainResult(matrices, residual, lookup, selected)


def solve_cell_full(
    internal_matrix: sparse.spmatrix,
) -> npt.NDArray[np.float64]:
    """Solve one cell's full reduction using the legacy exit calculation.

    Parameters
    ----------
    internal_matrix : scipy.sparse.spmatrix
        Square matrix Q of this cell's internal transition coefficients,
        including all global states. Supply COO format to preserve the legacy
        row-summation order.

    Returns
    -------
    numpy.ndarray
        Activation-by-exit-source matrix for this cell. Exit coefficients are
        calculated as one minus each row sum of Q.
    """
    exit_rates = 1 - np.asarray(internal_matrix.sum(axis=1)).ravel()
    system = (
        sparse.eye(internal_matrix.shape[0], format="csc") - internal_matrix
    ).tocsc()
    return splu(system, permc_spec="MMD_AT_PLUS_A").solve(np.diag(exit_rates))


def solve_cell_partial(
    internal_matrix: sparse.spmatrix,
    exit_rates: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Solve one cell's partial reduction from a specified subset of probabilities.

    Parameters
    ----------
    internal_matrix : scipy.sparse.spmatrix
        Square matrix Q of this cell's eliminated internal coefficients,
        in compact selected-state order.
    exit_rates : numpy.ndarray
        One-dimensional sums of this cell's retained channel coefficients,
        in the same selected-state order.

    Returns
    -------
    numpy.ndarray
        Activation-by-exit-source matrix B satisfying
        ``(I-Q) B = diag(exit_rates)`` for this cell.
    """
    system = (
        sparse.eye(internal_matrix.shape[0], format="csc") - internal_matrix
    ).tocsc()
    return splu(system, permc_spec="MMD_AT_PLUS_A").solve(np.diag(exit_rates))
