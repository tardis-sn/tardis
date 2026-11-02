"""Compiled-path verification of independent hybrid macro-atom visits."""

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest
from numba import njit

from tardis.opacities.continuum.continuum_state_numba import (
    ContinuumOpacityStateNumba,
)
from tardis.opacities.macro_atom.absorbing_markov_chain import (
    create_absorbing_probs,
)
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.montecarlo.macro_atom import macro_atom_interaction_iip


@njit
def sample_emissions(
    opacity: OpacityStateNumba,
    continuum: ContinuumOpacityStateNumba,
    activation: int,
    count: int,
    seed: int,
) -> npt.NDArray[np.int64]:
    """Draw independent packets through the production compiled interaction."""
    np.random.seed(seed)
    outcomes = np.zeros(3, dtype=np.int64)
    for _ in range(count):
        channel, _ = macro_atom_interaction_iip(
            activation, 0, opacity, continuum
        )
        outcomes[channel] += 1
    return outcomes


@pytest.fixture(
    params=[
        ("full", False, False),
        ("partial", False, False),
        ("partial", True, False),
        ("partial", False, True),
        ("partial", True, True),
        ("explicit", False, False),
    ]
)
def two_visit_state(request: pytest.FixtureRequest) -> tuple:
    """H visits k/i, with thermal emission and re-entry through unselected He."""
    mode, include_thermal_pool, include_ionization_pool = request.param
    metadata = pd.DataFrame(
        {
            "source_level_idx": [0, 0, 1, 2, 2, 3, 3, 4, 5, 6],
            "destination_level_idx": [1, -1, 2, 5, -1, 4, -1, -1, 6, 3],
            "transition_type": [1, -1, 11, 22, -21, 1, -7, -1, 3, 1],
            "source": [(1, 0, 0)] * 2
            + [(1, 0, 1), ("k", -99, -99), ("k", -99, -99)]
            + [(1, 1, 0)] * 2
            + [(1, 1, 1), ("i", -99, -99), (2, 0, 0)],
        }
    )
    result = create_absorbing_probs(
        pd.DataFrame([0.5, 0.5, 1.0, 0.75, 0.25, 0.5, 0.5, 1.0, 1.0, 1.0]),
        metadata,
        mode=mode,
        elements=["H"] if mode == "partial" else [],
        include_thermal_pool=include_thermal_pool,
        include_ionization_pool=include_ionization_pool,
    )
    residual = result.residual_probabilities.to_numpy(copy=True)
    totals = np.bincount(metadata.source_level_idx, weights=residual[:, 0])
    np.divide(
        residual[:, 0],
        totals[metadata.source_level_idx],
        out=residual[:, 0],
        where=totals[metadata.source_level_idx] != 0,
    )
    opacity = OpacityStateNumba(
        np.ones(1),
        np.ones(1),
        np.ones(1),
        np.ones((1, 1)),
        residual,
        np.array([0], dtype=np.int64),
        np.array([0, 2, 3, 5, 7, 8, 9, 10], dtype=np.int64),
        metadata.transition_type.to_numpy(dtype=np.int64, copy=True),
        metadata.destination_level_idx.to_numpy(dtype=np.int64),
        np.array([0, 0, 0, 0, 2, 0, 1, 2, 0, 0], dtype=np.int64),
    )
    continuum = ContinuumOpacityStateNumba(
        np.ones(1),
        np.ones((1, 1)),
        np.ones(1),
        np.ones(1),
        np.array([0, 1], dtype=np.int64),
        np.ones((1, 1)),
        np.ones(1),
        np.ones(1),
        np.ones(1),
        np.ones((1, 1)),
        np.array([5], dtype=np.int64),
        2,
        result.matrix,
        result.macro_level_to_reduced_idx,
        result.reduced_to_macro_level_idx,
    )
    return opacity, continuum


def test_compiled_reentry_emission_distribution(two_visit_state: tuple) -> None:
    """Exact tree law is 1/2, 3/16, 5/16 with thermal emission and re-entry.

    Tolerance is six binomial standard deviations per channel (60000 packets),
    giving a conservative familywise stochastic bound across all modes and pool selections.
    Reusing the first reduced random number suppresses one channel entirely.
    """
    opacity, continuum = two_visit_state
    count = 60000
    expected = np.array([0.5, 0.1875, 0.3125])
    actual = sample_emissions(opacity, continuum, 0, count, 1741) / count
    assert np.all(
        np.abs(actual - expected)
        <= 6 * np.sqrt(expected * (1 - expected) / count)
    ), actual
    # A k activation can emit directly; i goes through He before re-entering H.
    for activation, pool_expected in (
        (2, np.array([0.0, 0.375, 0.625])),
        (5, np.array([0.0, 0.5, 0.5])),
    ):
        pool_outcomes = (
            sample_emissions(opacity, continuum, activation, 20000, 972) / 20000
        )
        np.testing.assert_allclose(
            pool_outcomes,
            pool_expected,
            atol=6 * np.sqrt(0.25 / 20000),
            rtol=0,
        )
    assert macro_atom_interaction_iip.nopython_signatures
