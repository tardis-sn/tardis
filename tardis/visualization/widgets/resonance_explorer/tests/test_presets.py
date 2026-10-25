"""Checks for the manuscript's physical case families."""

from math import sqrt

import numpy.testing as npt
import pytest

from tardis.visualization.widgets.resonance_explorer.analysis import (
    analyze_case,
)
from tardis.visualization.widgets.resonance_explorer.presets import PRESETS


def test_presets_cover_distinct_geometry_and_solver_families() -> None:
    """Claim: each named single-shell family has a selectable exact input."""
    required = {
        "linked three crossings",
        "outbound one crossing",
        "closest approach launch",
        "radial outward",
        "inner direct hit",
        "inner grazing hit",
        "radial inward",
        "inner miss flyby",
        "inner tangent launch code limit",
        "static matched",
        "static missed",
        "uniform speed analytic",
        "uniform radial non-isolated",
        "uniform radial outward non-isolated",
        "radius proportional",
        "positive offset linear",
        "within-shell stagnation",
        "no turn radius",
        "turn outside shell",
        "stationary inflection flat center",
        "one forward turn two crossings",
        "target above range",
        "extremum touch",
        "scatter truncation",
        "launch endpoint excluded",
        "scatter before resonance",
        "scatter at resonance tie",
        "outer boundary target tie",
        "target below outbound range",
        "line center",
        "mirrored flow",
    }
    assert required <= PRESETS.keys()


def test_linked_preset_preserves_exact_turning_geometry() -> None:
    """Claim: the linked input gives the paper's ±0.9 gradient zeros."""
    inputs = PRESETS["linked three crossings"].inputs
    npt.assert_allclose(
        (inputs.r_in, inputs.r_out, inputs.b, inputs.w),
        (1.0, 3.0, 1.2, 0.3),
        atol=0,
        rtol=0,
    )
    npt.assert_allclose(inputs.x_0, -sqrt(7.56), atol=1e-15)
    slope = (inputs.v_out - inputs.v_in) / (inputs.r_out - inputs.r_in)
    intercept = inputs.v_in - slope * inputs.r_in
    npt.assert_allclose((slope, intercept), (1.0, -2.34375), atol=1e-15)
    npt.assert_allclose(
        -intercept * inputs.b**2 / slope,
        (inputs.b**2 + 0.9**2) ** 1.5,
        atol=1e-14,
    )


def test_scatter_preset_ends_at_stated_coordinate() -> None:
    """Claim: an assumed scatter clips the paper's later crossing."""
    inputs = PRESETS["scatter truncation"].inputs
    npt.assert_allclose(inputs.x_0 + inputs.s_stop, 1.6, atol=0)
    npt.assert_allclose(inputs.w, 0.275, atol=0)


@pytest.mark.parametrize("name", PRESETS)
def test_preset_matches_paper_topology(name: str) -> None:
    """Claim: named cases reproduce their independently cataloged topology."""
    preset = PRESETS[name]
    snapshot = analyze_case(preset.inputs)
    assert snapshot.ray_family == preset.expected_ray_family, name
    assert snapshot.shape_family == preset.expected_shape_family, name
    assert len(snapshot.turns) == preset.expected_turn_count, name
    if preset.expected_crossing_count is not None:
        assert len(snapshot.crossings) == preset.expected_crossing_count, name


def test_boundary_tie_preset_exposes_transport_ownership() -> None:
    """Claim: an outer-edge root is geometric but loses to boundary event."""
    snapshot = analyze_case(PRESETS["outer boundary target tie"].inputs)
    assert len(snapshot.crossings) == 1
    assert snapshot.nearest_resonance is None
    assert any(
        "boundary tie" in candidate.flags for candidate in snapshot.candidates
    )


def test_scatter_tie_preset_keeps_equal_distance_line_eligible() -> None:
    """Claim: strict scatter comparison preserves an equal-distance line."""
    snapshot = analyze_case(PRESETS["scatter at resonance tie"].inputs)
    assert any(
        candidate.event_eligible and candidate.s == snapshot.active_end
        for candidate in snapshot.candidates
    )
