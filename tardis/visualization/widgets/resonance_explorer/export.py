"""Export one immutable explorer analysis as figure or strict JSON bytes."""

from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import io
import json
import math
import subprocess
import sys
from dataclasses import fields, is_dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from .analysis import CaseAnalysis
from .diagnostics import diagnose_solver
from .evidence import build_evidence_record
from .fingerprint import LENGTH_UNIT, VELOCITY_UNIT, input_fingerprint
from .render import render_case

if TYPE_CHECKING:
    from matplotlib.figure import Figure


SCHEMA_VERSION = "1.3"


def json_safe(value: object) -> object:
    """Convert a snapshot field to values accepted by strict JSON."""
    if is_dataclass(value):
        safe_value = {
            field.name: json_safe(getattr(value, field.name))
            for field in fields(value)
        }
    elif isinstance(value, dict):
        safe_value = {str(key): json_safe(item) for key, item in value.items()}
    elif isinstance(value, (tuple, list, np.ndarray)):
        safe_value = [json_safe(item) for item in value]
    elif isinstance(value, np.generic):
        safe_value = json_safe(value.item())
    elif isinstance(value, float):
        safe_value = value if math.isfinite(value) else None
    elif isinstance(value, Path):
        safe_value = str(value)
    else:
        safe_value = value
    return safe_value


def _git_output(repository: Path, *arguments: str) -> str | None:
    """Read provenance from one checkout without changing it."""
    try:
        process = subprocess.run(
            ["git", "-C", str(repository), *arguments],
            capture_output=True,
            check=False,
            text=True,
            timeout=5,
        )
    except OSError, subprocess.TimeoutExpired:
        return None
    return process.stdout.strip() if process.returncode == 0 else None


def _package_version(name: str) -> str | None:
    """Read an installed version if this environment provides one."""
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def _provenance() -> dict[str, object]:
    """Record package and source state needed to reproduce an export."""
    explorer_sources = {}
    for source_path in sorted(Path(__file__).resolve().parent.glob("*.py")):
        source_bytes = source_path.read_bytes()
        explorer_sources[source_path.name] = {
            "sha256": hashlib.sha256(source_bytes).hexdigest(),
            "source_utf8": source_bytes.decode("utf-8"),
        }
    repository_root = Path(__file__).resolve().parents[3]
    solver_module = importlib.import_module(
        "tardis.transport.montecarlo.nonhomologous_grid"
    )
    solver_path = Path(solver_module.__file__).resolve()
    tardis_root = solver_path.parents[3]
    repository_status = _git_output(repository_root, "status", "--porcelain")
    return {
        "python_version": sys.version.split()[0],
        "package_versions": {
            name: _package_version(name)
            for name in ("tardis", "panel", "matplotlib", "numpy")
        },
        "tardis_module_path": str(solver_path),
        "tardis_source_commit": _git_output(tardis_root, "rev-parse", "HEAD"),
        "explorer_source_commit": _git_output(repository_root, "rev-parse", "HEAD"),
        "explorer_dirty": bool(repository_status) if repository_status is not None else None,
        "explorer_sources": explorer_sources,
    }


def build_export_record(analysis: CaseAnalysis) -> dict[str, object]:
    """Build a replayable record from one immutable analysis snapshot.

    Parameters
    ----------
    analysis : CaseAnalysis
        Snapshot used by the current figure and candidate table.

    Returns
    -------
    dict[str, object]
        JSON-safe scientific state and source provenance.
    """
    raw_slots = [
        {
            "slot_idx": slot_idx,
            "scaled_root": float(root) if math.isfinite(root) else None,
            "status": "finite" if math.isfinite(root) else "unavailable",
        }
        for slot_idx, root in enumerate(analysis.raw_roots)
    ]
    structured_parameters = getattr(analysis, "structured_parameters", ())
    solver = {
        "route": analysis.solver_route,
        "structured_parameters": json_safe(dict(structured_parameters)),
        "distance_tolerance": getattr(analysis, "distance_tolerance", None),
        "raw_slots": raw_slots,
        "candidates": json_safe(analysis.candidates),
    }
    if getattr(analysis, "solver_provenance", None) is not None:
        solver["reconstructed_diagnostic"] = json_safe(
            diagnose_solver(analysis)
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "input_fingerprint": input_fingerprint(analysis.inputs),
        "units": {"length": LENGTH_UNIT, "velocity": VELOCITY_UNIT},
        "inputs": json_safe(analysis.inputs),
        "case": {
            "ray_family": analysis.ray_family,
            "shape_family": analysis.shape_family,
            "result_family": analysis.result_family,
        },
        "solver": json_safe(solver),
        "analysis": json_safe(analysis),
        "evidence": json_safe(build_evidence_record(analysis))
        if isinstance(analysis, CaseAnalysis)
        else None,
        "provenance": _provenance(),
    }


def export_json_bytes(analysis: CaseAnalysis) -> bytes:
    """Serialize an analysis snapshot to strict UTF-8 JSON bytes.

    Parameters
    ----------
    analysis : CaseAnalysis
        Current committed explorer snapshot.

    Returns
    -------
    bytes
        Indented JSON with nonfinite values represented by ``null``.
    """
    payload = build_export_record(analysis)
    return (
        json.dumps(payload, allow_nan=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")


def export_figure_bytes(
    analysis: CaseAnalysis,
    format_name: str,
    selected_slot_index: int | None = None,
) -> bytes:
    """Render PDF or PNG bytes from the same committed analysis snapshot.

    Parameters
    ----------
    analysis : CaseAnalysis
        Current committed explorer snapshot.
    format_name : str
        ``"pdf"`` or ``"png"``.
    selected_slot_index : int or None, optional
        Raw solver slot to emphasize on the figure.

    Returns
    -------
    bytes
        Encoded figure data.
    """
    if format_name not in {"pdf", "png"}:
        raise ValueError("format_name must be 'pdf' or 'png'")
    figure: Figure = render_case(analysis, selected_slot_index)
    output = io.BytesIO()
    try:
        figure.savefig(output, format=format_name, dpi=180, bbox_inches="tight")
    finally:
        figure.clear()
    return output.getvalue()
