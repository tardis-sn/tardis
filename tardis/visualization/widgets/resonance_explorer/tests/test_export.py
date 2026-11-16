"""Checks for replayable and strict snapshot export."""

# ruff: noqa: TID252 - Package-local tests use relative imports under pytest.

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

from ..analysis import analyze_case
from ..export import (
    export_json_bytes,
)
from ..fingerprint import (
    input_fingerprint,
)
from ..model import ExplorerInputs


@dataclass(frozen=True)
class ExampleInputs:
    r_in: float = 1.0
    r_out: float = 3.0
    v_in: float = -1.34375
    v_out: float = 0.65625
    b: float = 1.2
    x_0: float = -math.sqrt(9.0 - 1.2**2)
    w: float = 0.3
    s_stop: float | None = None
    stop_cause: str = "assumed electron scatter"


@dataclass(frozen=True)
class ExampleAnalysis:
    inputs: ExampleInputs
    m: float = 1.0
    q: float = -2.34375
    solver_route: str = "general quartic"
    raw_roots: tuple[float, ...] = (math.nan, -1.0, 2.0, math.nan)
    structured_parameters: tuple[tuple[str, float], ...] = (("N", 0.15),)
    distance_tolerance: float = 1.0e-12
    candidates: tuple[object, ...] = ()
    ray_family: str = "inner-miss flyby"
    shape_family: str = "paired projected turns"
    result_family: str = "three forward crossings"


def test_input_fingerprint_tracks_exact_inputs() -> None:
    """Changing a scientific control must change the exported state identity."""
    first = ExampleInputs()
    changed = ExampleInputs(w=0.7)
    assert input_fingerprint(first) == input_fingerprint(first)
    assert len(input_fingerprint(first)) == 16
    assert input_fingerprint(first) != input_fingerprint(changed)


def test_json_export_retains_all_raw_slots_without_nonfinite_literals() -> None:
    """NaN production slots remain addressable with an explicit status."""
    analysis = ExampleAnalysis(inputs=ExampleInputs())
    encoded = export_json_bytes(analysis)
    decoded = json.loads(
        encoded,
        parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)),
    )

    assert decoded["schema_version"] == "1.3"
    assert decoded["input_fingerprint"] == input_fingerprint(analysis.inputs)
    assert decoded["inputs"]["w"] == 0.3
    assert [slot["slot_idx"] for slot in decoded["solver"]["raw_slots"]] == [
        0,
        1,
        2,
        3,
    ]
    assert decoded["solver"]["raw_slots"][0] == {
        "slot_idx": 0,
        "scaled_root": None,
        "status": "unavailable",
    }
    assert decoded["solver"]["raw_slots"][1]["scaled_root"] == -1.0
    assert decoded["solver"]["route"] == "general quartic"
    assert decoded["case"]["ray_family"] == "inner-miss flyby"
    assert decoded["provenance"]["tardis_source_commit"]


def test_real_export_embeds_independent_evidence_and_source_revision() -> None:
    """A replay record carries both original-equation roots and source hashes."""
    inputs = ExplorerInputs(1, 3, -1.34375, 0.65625, 1.2, -math.sqrt(7.56), 0.3)
    decoded = json.loads(export_json_bytes(analyze_case(inputs)))
    evidence = decoded["evidence"]

    assert len(evidence["independent_roots"]) == 3
    assert len(evidence["production_slots"]) == 4
    assert evidence["route"] == decoded["solver"]["route"]
    assert (
        evidence["source_revision"]["tardis_commit"]
        == decoded["provenance"]["tardis_source_commit"]
    )
    assert (
        evidence["source_revision"]["explorer_sha256"]["evidence.py"]
        == (decoded["provenance"]["explorer_sources"]["evidence.py"]["sha256"])
    )


def test_json_export_embeds_reconstructable_runtime_source(
    tmp_path: Path,
) -> None:
    """Claim: an export can restore the exact explorer runtime source.

    Verification: compare restored bytes and SHA-256 with every top-level
    Python source file used by the explorer. This verifies source provenance
    wiring; scientific behavior is checked separately.
    """
    package_directory = Path(__file__).resolve().parents[1]
    original_paths = sorted(package_directory.glob("*.py"))
    decoded = json.loads(export_json_bytes(ExampleAnalysis(ExampleInputs())))
    sources = decoded["provenance"]["explorer_sources"]

    assert set(sources) == {path.name for path in original_paths}
    for original_path in original_paths:
        saved = sources[original_path.name]
        restored_path = tmp_path / original_path.name
        restored_path.write_bytes(saved["source_utf8"].encode("utf-8"))
        original_bytes = original_path.read_bytes()
        restored_bytes = restored_path.read_bytes()
        assert restored_bytes == original_bytes
        assert saved["sha256"] == hashlib.sha256(restored_bytes).hexdigest()
