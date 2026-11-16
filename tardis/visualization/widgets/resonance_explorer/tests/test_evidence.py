"""Independent original-equation and packet-distance evidence checks."""

import math

import numpy as np
import pytest

from .. import evidence as evidence_module  # noqa: TID252
from ..analysis import analyze_case  # noqa: TID252
from ..evidence import build_evidence_record  # noqa: TID252
from ..model import ExplorerInputs  # noqa: TID252


def _case(m: float, q: float, w: float, x_0: float) -> ExplorerInputs:
    """Construct a common-shell input with explicit m and q."""
    return ExplorerInputs(1, 3, m + q, 3 * m + q, 1.2, x_0, w)


@pytest.mark.parametrize(
    ("m", "q", "w", "x_0", "expected_x"),
    [
        (1, 0, 0.3, -math.sqrt(7.56), [0.3]),
        (
            1,
            -2.34375,
            0.3,
            -math.sqrt(7.56),
            [-1.55588727, -0.341472541, 2.395531338],
        ),
        (1, -2.34375, 0.3, -0.5, [-0.341472541, 2.395531338]),
        (1, -2.34375, 0.3, 1.2, [2.395531338]),
        (1, -3.5, 0.6, -math.sqrt(7.56), [-2.571701094, -0.331209962]),
        (-1, 0, 0.6, -math.sqrt(7.56), [-0.6]),
    ],
)
def test_guided_independent_roots_match_analytic_or_fixed_brackets(
    m: float, q: float, w: float, x_0: float, expected_x: list[float]
) -> None:
    """Claim: guided crossings satisfy independently specified unsquared roots.

    Regime: common full or clipped flyby, including wholly collapsing shell.
    Verification: homology has analytic x=w/m; nonlinear expected values are
    independently bracketed coordinates printed to 9 decimal places in the
    brief or prior manuscript audit. The 2e-9 tolerance covers that printing.
    """
    evidence = build_evidence_record(analyze_case(_case(m, q, w, x_0)))
    roots = evidence["independent_roots"]
    np.testing.assert_allclose(
        [root["x"] for root in roots], expected_x, atol=2e-9, rtol=0
    )
    assert all(abs(root["residual"]) < 1e-12 for root in roots)
    assert evidence["agreement"]["status"] == "matched"
    assert len(evidence["production_slots"]) == 4


@pytest.mark.parametrize(
    ("m", "q", "w", "expected_distance"),
    [
        (1, 0, 0.3, 0.3 + math.sqrt(7.56)),
        (1, -2.34375, 0.3, -1.55588727 + math.sqrt(7.56)),
        (1, -3.5, 0.6, -2.571701094 + math.sqrt(7.56)),
        (-1, 0, 0.6, -0.6 + math.sqrt(7.56)),
    ],
)
def test_packet_nearest_agrees_with_analytic_or_bracketed_distance(
    m: float, q: float, w: float, expected_distance: float
) -> None:
    """Claim: ordinary, repeated-crossing, and collapse packet distances agree.

    Verification: the expected distance is a separately known unsquared root
    minus the fixed launch coordinate. Scale mapping is 1e14 cm per length
    unit and 1e8 cm/s per velocity unit. The 2e-9 length tolerance includes
    the 9-decimal bracket print precision and frequency-to-target rounding.
    """
    evidence = build_evidence_record(
        analyze_case(_case(m, q, w, -math.sqrt(7.56)))
    )
    packet = evidence["packet_nearest"]
    assert packet["status"] == "compared"
    assert packet["comparison_valid"]
    np.testing.assert_allclose(
        packet["packet_distance"], expected_distance, atol=2e-9, rtol=0
    )
    assert packet["agreement_within_1e-8"]


def test_linked_fourth_slot_has_wrong_unsquared_sign() -> None:
    """Claim: four quartic slots do not imply four physical crossings.

    Verification: the linked case's three independently bracketed roots
    coexist with one production slot rejected by the original sign relation.
    """
    evidence = build_evidence_record(
        analyze_case(_case(1, -2.34375, 0.3, -math.sqrt(7.56)))
    )
    assert len(evidence["independent_roots"]) == 3
    assert (
        sum(
            "squared-sign artifact" in slot["flags"]
            for slot in evidence["production_slots"]
        )
        == 1
    )


def test_exact_contact_mapping_reports_frequency_rounding_limit() -> None:
    """Claim: a mapped exact tangent is not an exact packet parity check.

    The 1e15 Hz line frequency rounds the recovered velocity target 4.8e-14
    above the analytic crest. The packet therefore selects the outbound simple
    root, while the illustrative exact target has an earlier repeated contact.
    The visible mismatch must not be presented as a solver parity failure.
    """
    evidence = build_evidence_record(
        analyze_case(_case(1, -2.34375, 0.50625, -math.sqrt(7.56)))
    )
    roots = evidence["independent_roots"]
    np.testing.assert_allclose(
        [root["x"] for root in roots], [-0.9, 2.639908811], atol=2e-9
    )
    packet = evidence["packet_nearest"]
    assert packet["status"] == "contact sensitive to frequency rounding"
    assert not packet["comparison_valid"]
    assert packet["recovered_target"] > 0.50625
    assert packet["difference"] > 3.5


def test_evidence_source_hashes_survive_unavailable_git(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Source hashes remain available when Git cannot report a commit."""

    def unavailable_git(*_args: object, **_kwargs: object) -> None:
        raise OSError("git unavailable")

    monkeypatch.setattr(evidence_module.subprocess, "run", unavailable_git)
    revision = evidence_module._source_revision()
    assert revision["explorer_commit"] == revision["tardis_commit"]
    assert revision["tardis_commit"] is None
    assert revision["explorer_sha256"]["evidence.py"]
    assert revision["transport_sha256"]
