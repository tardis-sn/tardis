"""Verify direct continuum traversal across the split opacity interface."""

from types import SimpleNamespace

import numpy as np
import numpy.testing as npt
import pandas as pd

from tardis.opacities.continuum.continuum_state import ContinuumOpacityState
from tardis.opacities.opacity_state import OpacityState
from tardis.transport.montecarlo.macro_atom import MacroAtomTransitionType
from tardis.transport.montecarlo.modes.iip.interaction_event_callers import (
    macro_atom_event,
)
from tardis.transport.montecarlo.packets.radiative_packet import RPacket


def test_opacity_conversion_preserves_optical_depth_coefficient() -> None:
    """Preserve the supplied line coefficient through Numba conversion."""
    opacity = OpacityState(
        pd.Series([1.0]),
        np.array([8000.0]),
        pd.Series([2e14]),
        pd.DataFrame([[1.0]]),
        pd.DataFrame([[1.0]]),
        sobolev_optical_depth_coefficient=pd.DataFrame([[7.0]]),
    )
    actual = opacity.to_numba(None, "scatter")
    npt.assert_array_equal(actual.sobolev_optical_depth_coefficient, [[7.0]])


def test_direct_continuum_event_follows_internal_transition() -> None:
    """An internal jump followed by certain line emission has a known frequency.

    Regime: two macro levels in one static shell, with deterministic 0->1
    internal transition and 1->line deactivation. Empty reduction maps retain
    explicit traversal. Verification uses unit branch probabilities
    and the prescribed emitted frequency 2e14 Hz, not a saved implementation
    output. Relative tolerance 1e-14 allows floating-point frame transforms.
    This also compiles the actual IIP event handler with split Numba states.
    """
    metadata = pd.DataFrame(
        {
            "transition_type": [
                MacroAtomTransitionType.INTERNAL_UP,
                MacroAtomTransitionType.BB_EMISSION,
            ],
            "destination_level_idx": [1, 0],
            "transition_line_idx": [0, 0],
        },
        dtype=np.int64,
    )
    macro_state = SimpleNamespace(
        transition_probabilities=pd.DataFrame([[1.0], [1.0]]),
        normalized_deactivating_probs=pd.DataFrame([[1.0], [1.0]]),
        line2macro_level_upper=pd.Series([0], dtype=np.int64),
        macro_block_edge_index=np.array([0, 1, 2], dtype=np.int64),
        transition_metadata=metadata,
        absorbing_probability_matrix=np.empty((1, 0, 0)),
        macro_level_to_reduced_idx=np.full(3, -1, dtype=np.int64),
        reduced_to_macro_level_idx=np.empty(0, dtype=np.int64),
        photo_ion_block_idx=-1,
        k_packet_idx=2,
    )
    line_opacity = OpacityState(
        pd.Series([1.0]),
        np.array([8000.0]),
        pd.Series([2e14]),
        pd.DataFrame([[1.0]]),
        pd.DataFrame([[1.0]]),
    )
    level_index = pd.MultiIndex.from_tuples([(2, 0, 0)])
    continuum = ContinuumOpacityState(
        pd.Series([1e15], index=level_index),
        pd.Series([0], index=level_index),
        pd.DataFrame([[1.0]], index=level_index),
        pd.DataFrame(
            {"nu": [1e15, 2e15], "x_sect": [0.0, 0.0]},
            index=level_index.repeat(2),
        ),
        pd.DataFrame([[0.0], [0.0]], index=level_index.repeat(2)),
        np.array([0.0]),
        pd.DataFrame([[0.0], [1.0]], index=level_index.repeat(2)),
        pd.DataFrame({"destination_level_idx": [0]}, index=level_index),
    )
    lines_numba = line_opacity.to_numba(
        macro_state, "macroatom", continuum_processes_enabled=True
    )
    continuum_numba = continuum.to_numba(np.array([8000.0]), macro_state)
    assert continuum_numba.absorbing_markov_probabilities.shape == (1, 0, 0)
    packet = RPacket(0.0, 0.0, 1e14, 1.0, 123)
    macro_atom_event(0, packet, 1e6, lines_numba, False, continuum_numba)
    assert macro_atom_event.nopython_signatures
    npt.assert_allclose(packet.nu, 2e14, rtol=1e-14, atol=0.0)
    npt.assert_allclose(packet.energy, 1.0, rtol=1e-14, atol=0.0)
