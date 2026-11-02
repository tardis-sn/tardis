"""Calculate first-exit probabilities for reduced macro-atom chains.

This module eliminates selected internal transitions independently in each cell
and computes the source-state probabilities for the first retained transition.
It does not calculate expected step counts or visit counts, though the formalism
does allow for that. It would require an additional potentially expensive linear
algebra solve so it has been omitted.

References
----------
    Absorbing Markov chain theory: https://en.wikipedia.org/wiki/Absorbing_Markov_chain
"""

from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal, localcontext
from functools import partial
from typing import NamedTuple

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import splu

from tardis.util.base import element_symbol2atomic_number

# LU/CDF arithmetic tolerance, not a tolerance for discarding physical channels.
PROBABILITY_TOLERANCE = 1e-10


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
    *,
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
    transition_probabilities : pandas.DataFrame
        Normalized original channel coefficients; columns are cells.
    metadata : pandas.DataFrame
        Transition source tuples, global source/destination indices, and types.
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
        Cell-by-activation-by-exit-source matrix, residual channel probabilities,
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
        selected_rows = metadata.source.map(
            lambda level: level[0] in atomic_numbers
        ).to_numpy()
        represented = {level[0] for level in metadata.source}
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

    cell_solver = partial(
        solve_cell,
        transition_probabilities=transition_probabilities,
        transition_index=metadata.index,
        sources=sources,
        destinations=destinations,
        internal=internal,
        global_size=global_size,
        selected=selected,
        eliminated=eliminated,
        residual=residual,
        rows=rows,
        columns=columns,
        mode=mode,
    )

    if max_workers > 1 and len(matrices) > 1:
        with ThreadPoolExecutor(
            max_workers=min(max_workers, len(matrices))
        ) as executor:
            for cell, solved in executor.map(cell_solver, range(len(matrices))):
                matrices[cell] = solved
    else:
        for cell in range(len(matrices)):
            _, matrices[cell] = cell_solver(cell)
    return AbsorbingChainResult(matrices, residual, lookup, selected)


def solve_cell(
    cell: int,
    *,
    transition_probabilities: pd.DataFrame,
    transition_index: pd.Index,
    sources: npt.NDArray[np.int64],
    destinations: npt.NDArray[np.int64],
    internal: npt.NDArray[np.bool_],
    global_size: int,
    selected: npt.NDArray[np.int64],
    eliminated: npt.NDArray[np.bool_],
    residual: pd.DataFrame,
    rows: npt.NDArray[np.int64],
    columns: npt.NDArray[np.int64],
    mode: str,
) -> tuple[int, npt.NDArray[np.float64]]:
    """Check one cell's physical exits and solve its selected first-exit law.

    Parameters
    ----------
    cell : int
        Positional cell column to solve.
    transition_probabilities : pandas.DataFrame
        Normalized original channel coefficients for all cells.
    transition_index : pandas.Index
        Original transition labels for diagnostics.
    sources, destinations : numpy.ndarray
        Global source and destination indices for every original channel.
    internal : numpy.ndarray
        Mask identifying original internal transitions.
    global_size : int
        Number of global level indices, including gaps.
    selected : numpy.ndarray
        Global states in compact matrix order.
    eliminated : numpy.ndarray
        Mask of selected-to-selected internal channels to eliminate.
    residual : pandas.DataFrame
        Original channels with eliminated coefficients set to zero.
    rows, columns : numpy.ndarray
        Compact source and destination indices for eliminated channels.
    mode : {'full', 'partial', 'explicit'}
        Reduction mode, including legacy full-mode numerical behavior.

    Returns
    -------
    tuple[int, numpy.ndarray]
        Cell position and its first-exit matrix; explicit mode returns an
        empty matrix after checking the original channels.

    Notes
    -----
    Shared input tables and index arrays are read without modification so
    cells can be solved independently in parallel.
    """
    values = transition_probabilities.iloc[:, cell].to_numpy()
    totals = np.bincount(sources, weights=values, minlength=global_size)
    active = totals > 0
    if not np.all(np.isfinite(values)):
        raise ValueError(
            f"cell {cell}: nonfinite channel coefficients at transitions {transition_index[~np.isfinite(values)].tolist()}"
        )
    if not np.allclose(totals[active], 1, rtol=0, atol=PROBABILITY_TOLERANCE):
        raise ValueError(
            f"cell {cell}: original outgoing probabilities are not normalized"
        )
    # Reverse reachability from positive physical deactivation channels detects
    # closed communicating classes even in explicit mode, before MC can trap.
    reachable = (
        np.bincount(
            sources[~internal],
            weights=values[~internal],
            minlength=global_size,
        )
        > 0
    )
    positive_internal = internal & (values > 0)
    graph = sparse.csr_matrix(
        (
            np.ones(positive_internal.sum()),
            (destinations[positive_internal], sources[positive_internal]),
        ),
        shape=(global_size, global_size),
    )
    pending = np.flatnonzero(reachable).tolist()
    while pending:
        destination = pending.pop()
        for source in graph.indices[
            graph.indptr[destination] : graph.indptr[destination + 1]
        ]:
            if not reachable[source]:
                reachable[source] = True
                pending.append(source)
    trapped = np.flatnonzero(active & ~reachable)
    if len(trapped):
        raise ValueError(
            f"cell {cell}: closed or non-emitting communicating class reachable from states {trapped}"
        )
    if not len(selected):
        return cell, np.zeros((0, 0), dtype=np.float64)
    q_matrix = sparse.coo_matrix(
        (values[eliminated], (rows, columns)),
        shape=(len(selected), len(selected)),
    )
    exit_rates = np.bincount(
        sources,
        weights=residual.iloc[:, cell].to_numpy(),
        minlength=global_size,
    )[selected]
    q_totals = np.asarray(q_matrix.sum(axis=1)).ravel()
    if not np.allclose(
        q_totals + exit_rates,
        totals[selected],
        rtol=0,
        atol=PROBABILITY_TOLERANCE,
    ):
        raise ValueError(
            f"cell {cell}: eliminated and residual channels do not partition the original rates"
        )
    system = (sparse.eye(len(selected), format="csc") - q_matrix).tocsc()
    rhs = exit_rates
    if mode == "full":
        rhs = 1 - q_totals
    # For near-trapping self-loops, form the diagonal from physical exits
    # and off-diagonal jumps.
    near_trapping = q_matrix.diagonal() > 0.99
    if mode != "full" and np.any(near_trapping):
        off_diagonal = q_matrix.copy()
        off_diagonal.setdiag(0)
        stable_diagonal = (
            np.asarray(off_diagonal.sum(axis=1)).ravel() + exit_rates
        )
        diagonal = system.diagonal()
        diagonal[near_trapping] = stable_diagonal[near_trapping]
        system.setdiag(diagonal)
        rhs = rhs.copy()
        rhs[near_trapping] = exit_rates[near_trapping]
    try:
        solved = splu(system, permc_spec="MMD_AT_PLUS_A").solve(np.diag(rhs))
    except RuntimeError as error:
        if mode == "full":
            raise ValueError(
                f"cell {cell}: closed selected chain at states {selected}"
            ) from error
        solved = solve_near_trapping_chain(
            q_matrix, exit_rates, active[selected]
        )
    if mode != "full" and not np.allclose(
        solved[active[selected]].sum(axis=1),
        1,
        rtol=0,
        atol=32 * np.finfo(np.float64).eps,
    ):
        solved = solve_near_trapping_chain(
            q_matrix, exit_rates, active[selected]
        )
    residual_norm = np.max(np.abs(system @ solved - np.diag(rhs)))
    residual_scale = max(
        1.0,
        sparse.linalg.norm(system, np.inf) * np.linalg.norm(solved, np.inf),
    )
    if (
        not np.all(np.isfinite(solved))
        or residual_norm > PROBABILITY_TOLERANCE * residual_scale
        or not np.allclose(
            solved[active[selected]].sum(axis=1),
            1,
            rtol=0,
            atol=PROBABILITY_TOLERANCE,
        )
    ):
        raise ValueError(
            f"cell {cell}: invalid absorbing solve for states {selected}; residual {residual_norm}"
        )
    return cell, solved


def solve_near_trapping_chain(
    internal_matrix: sparse.spmatrix,
    exit_rates: npt.NDArray[np.float64],
    active: npt.NDArray[np.bool_],
) -> npt.NDArray[np.float64]:
    """Resolve exceptional cancellation in a transient sparse LU solve.

    Parameters
    ----------
    internal_matrix : scipy.sparse.spmatrix
        Selected internal probabilities Q.
    exit_rates : numpy.ndarray
        Actual retained channel sums, without subtracting from one.
    active : numpy.ndarray
        States having original outgoing probability.

    Returns
    -------
    numpy.ndarray
        First-exit matrix evaluated with at least 80-digit elimination, then rounded
        to float64. Precision also resolves the smallest positive exit with 40
        guard digits. Used only if LU violates probability normalization.

    Notes
    -----
    The conservation-form diagonal is the sum of off-diagonal internal jumps
    and actual exits. This preserves tiny exits in multi-state near-closed
    cycles, where even a long residence time must give total exit probability
    one. Inactive rows represent no physical activation and stay zero.
    """
    size = len(exit_rates)
    with localcontext() as context:
        positive_exit_rates = exit_rates[exit_rates > 0]
        context.prec = max(
            80, 40 - int(np.floor(np.log10(positive_exit_rates.min())))
        )
        zero = Decimal(0)
        one = Decimal(1)
        dense = internal_matrix.toarray()
        system = [
            [
                -Decimal(float(dense[row, column])) if row != column else zero
                for column in range(size)
            ]
            for row in range(size)
        ]
        rhs = [[zero for _ in range(size)] for _ in range(size)]
        for row in range(size):
            rate = Decimal(float(exit_rates[row]))
            system[row][row] = (
                sum(
                    (
                        -system[row][column]
                        for column in range(size)
                        if row != column
                    ),
                    rate,
                )
                if active[row]
                else one
            )
            rhs[row][row] = rate
        for pivot_idx in range(size):
            pivot_row = max(
                range(pivot_idx, size),
                key=lambda row: abs(system[row][pivot_idx]),
            )
            if not system[pivot_row][pivot_idx]:
                raise ValueError(
                    "Closed selected chain in high-precision solve"
                )
            system[pivot_idx], system[pivot_row] = (
                system[pivot_row],
                system[pivot_idx],
            )
            rhs[pivot_idx], rhs[pivot_row] = rhs[pivot_row], rhs[pivot_idx]
            pivot = system[pivot_idx][pivot_idx]
            for row in range(pivot_idx + 1, size):
                factor = system[row][pivot_idx] / pivot
                for column in range(pivot_idx + 1, size):
                    system[row][column] -= factor * system[pivot_idx][column]
                for column in range(size):
                    rhs[row][column] -= factor * rhs[pivot_idx][column]
        solved = [[zero for _ in range(size)] for _ in range(size)]
        for row in range(size - 1, -1, -1):
            for channel in range(size):
                solved[row][channel] = (
                    rhs[row][channel]
                    - sum(
                        (
                            system[row][column] * solved[column][channel]
                            for column in range(row + 1, size)
                        ),
                        zero,
                    )
                ) / system[row][row]
        return np.asarray(solved, dtype=np.float64)
