"""Orchestration checks for the local explorer interface."""

# ruff: noqa: TID252 - Package-local tests use relative imports.

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from importlib import import_module

import pytest
from bokeh.document import Document
from matplotlib.figure import Figure

from tardis.transport.montecarlo.configuration.constants import (
    CLOSE_LINE_THRESHOLD,
)

from ..app import ResonanceExplorerWidget
from ..lessons import LESSONS
from ..model import ExplorerInputs
from ..presets import (
    DEFAULT_PRESET,
    PRESETS,
)


@dataclass(frozen=True)
class ExampleAnalysis:
    inputs: ExplorerInputs
    m: float = 1.0
    q: float = -2.0
    r_0: float = 2.0
    mu_0: float = 0.0
    boundary_distance: float = 1.0
    boundary_kind: str = "outer"
    active_end: float = 1.0
    x_end: float = 1.0
    stop_reason: str = "outer boundary"
    turns: tuple[object, ...] = ()
    candidates: tuple[object, ...] = ()
    crossings: tuple[object, ...] = ()
    nearest_resonance: object = None
    solver_route: str = "general quartic"
    ray_family: str = "outbound shell leg"
    shape_family: str = "monotone projection"
    result_family: str = "target miss"
    radial_law: str = "offset linear"


def test_invalid_edit_clears_previous_snapshot_and_all_numerical_views() -> (
    None
):
    """A rejected input must not leave a stale plot or table beside the error."""
    original = ExplorerInputs(1.0, 3.0, -1.0, 1.0, 1.2, 1.5, 0.3)

    def analyze(inputs: ExplorerInputs) -> ExampleAnalysis:
        if inputs.r_out <= inputs.r_in:
            raise ValueError("outer radius must exceed inner radius")
        return ExampleAnalysis(inputs=inputs)

    application = ResonanceExplorerWidget(
        initial_inputs=original,
        analyzer=analyze,
        renderer=lambda _analysis, _selection: Figure(),
    )
    assert application.analysis is not None
    assert application.figure_pane.object is not None

    application.set_inputs(replace(original, r_out=0.5))

    assert application.analysis is None
    assert application.figure_pane.object is None
    assert len(application.candidate_table.value) == 0
    assert (
        "outer radius must exceed inner radius" in application.error_pane.object
    )
    assert application.pdf_download.disabled
    assert application.png_download.disabled
    assert application.json_download.disabled


def test_default_case_serializes_to_a_browser_document() -> None:
    """The advanced diagnostic must contain only values Panel can serialize."""
    application = ResonanceExplorerWidget()
    document = Document()

    application.layout.server_doc(document)

    assert document.roots
    assert application.diagnostic_pane.object["provenance"].startswith(
        "reconstructed"
    )


def test_empty_numeric_widget_clears_results_and_recovers_on_valid_edit() -> (
    None
):
    """An in-progress field edit must not retain the previous scientific case."""
    application = ResonanceExplorerWidget()
    assert application.analysis is not None
    assert application.figure_pane.object is not None
    assert len(application.candidate_table.value) > 0

    application.controls["w"].value = None

    assert application.analysis is None
    assert application.figure_pane.object is None
    assert (
        application.summary_pane.object
        == "No current analysis: correct the input."
    )
    assert len(application.candidate_table.value) == 0
    assert application.diagnostic_pane.object == {}
    assert application.error_pane.visible
    assert "Line target w" in application.error_pane.object
    assert application.preset_select.value is None
    assert "Custom / modified inputs" in application.preset_detail.object
    assert application.pdf_download.disabled
    assert application.png_download.disabled
    assert application.json_download.disabled

    application.controls["w"].value = 0.7

    assert application.analysis is not None
    assert application.analysis.inputs.w == 0.7
    assert application.figure_pane.object is not None
    assert len(application.candidate_table.value) > 0
    assert application.diagnostic_pane.object
    assert not application.error_pane.visible
    assert not application.pdf_download.disabled
    assert not application.png_download.disabled
    assert not application.json_download.disabled


def test_high_residual_slot_is_visible_without_claiming_a_verified_crossing() -> (
    None
):
    """Show production acceptance separately from unsquared verification."""
    inputs = ExplorerInputs(
        1,
        3,
        -1103.1408330241957,
        -1103.137596530609,
        1.0987592278970888,
        -2.7915458368278663,
        22.32587042642611,
    )
    application = ResonanceExplorerWidget(initial_inputs=inputs)

    frame = application.candidate_table.value
    accepted = frame.loc[frame["production accepted"]]
    assert len(accepted) == 1
    slot = accepted.iloc[0]
    assert slot["production event eligible"]
    assert not slot["unsquared verified"]
    assert not slot["verified event eligible"]
    assert abs(slot["residual"]) > slot["unsquared tolerance"]

    summary = application.summary_pane.object
    assert "distinct verified forward crossings: 0" in summary
    assert "Production accepted but unsquared unverified: 1" in summary
    assert "Production-selected line candidate: **slot 1" in summary
    assert "Nearest verified eligible resonance: **none**" in summary


@pytest.mark.parametrize(
    ("field_name", "edited_value"),
    [
        ("r_in", 0.9),
        ("r_out", 3.1),
        ("v_in", -1.2),
        ("v_out", 0.7),
        ("b", 1.1),
        ("x_0", -2.6),
        ("w", 0.7),
    ],
)
def test_each_numeric_edit_marks_custom_and_reset_restores_preset(
    field_name: str, edited_value: float
) -> None:
    """Every numeric knob changes the identity of the selected named case."""
    application = ResonanceExplorerWidget(
        analyzer=ExampleAnalysis,
        renderer=lambda _analysis, _selection: Figure(),
    )
    original = PRESETS[DEFAULT_PRESET].inputs

    application.controls[field_name].value = edited_value

    assert application.preset_select.value is None
    assert "Custom / modified inputs" in application.preset_detail.object
    assert getattr(application.analysis.inputs, field_name) == edited_value

    application._on_reset(None)

    assert application.preset_select.value == DEFAULT_PRESET
    assert application.analysis.inputs == original
    assert all(
        application.controls[name].value == getattr(original, name)
        for name in application.controls
    )
    assert "linked three crossings" in application.preset_detail.object


def test_stop_controls_mark_custom_and_reset_restores_optional_stop() -> None:
    """Enabling and editing an assumed scatter cannot retain preset identity."""
    application = ResonanceExplorerWidget(
        analyzer=ExampleAnalysis,
        renderer=lambda _analysis, _selection: Figure(),
    )

    application.stop_enabled.value = True
    assert application.preset_select.value is None
    assert application.analysis.inputs.s_stop == 0.5
    application.stop_distance.value = 0.4
    assert application.analysis.inputs.s_stop == 0.4
    assert "Custom / modified inputs" in application.preset_detail.object

    application._on_reset(None)
    assert application.preset_select.value == DEFAULT_PRESET
    assert not application.stop_enabled.value
    assert application.stop_distance.disabled
    assert application.analysis.inputs.s_stop is None


def test_reselecting_named_scatter_case_restores_edited_inputs_and_detail() -> (
    None
):
    """A named preset restores its optional stop and recomputed detail card."""
    application = ResonanceExplorerWidget(renderer=lambda _analysis, _selection: Figure())
    key = "scatter truncation"
    preset = PRESETS[key]

    application.preset_select.value = key
    assert application.analysis.inputs == preset.inputs
    assert application.stop_enabled.value
    assert application.stop_distance.value == preset.inputs.s_stop

    application.controls["w"].value = preset.inputs.w + 0.01
    assert application.preset_select.value is None
    assert "Custom / modified inputs" in application.preset_detail.object

    application.preset_select.value = key
    assert application.analysis.inputs == preset.inputs
    assert application.controls["w"].value == preset.inputs.w
    assert application.stop_distance.value == preset.inputs.s_stop
    assert application.stop_enabled.value
    assert preset.name in application.preset_detail.object
    assert f"s={repr(application.analysis.active_end)}" in (
        application.preset_detail.object
    )


def test_every_named_preset_detail_uses_its_computed_snapshot() -> None:
    """The selected card exposes exact inputs, derived values, and raw decisions."""
    application = ResonanceExplorerWidget(renderer=lambda _analysis, _selection: Figure())

    for key, preset in PRESETS.items():
        application.preset_select.value = key
        snapshot = application.analysis
        detail = application.preset_detail.object
        assert snapshot.inputs == preset.inputs, key
        assert preset.name in detail, key
        for name in ("r_in", "r_out", "b", "x_0", "w"):
            assert f"{name}={repr(getattr(preset.inputs, name))}" in detail, key
        for name in ("m", "q"):
            assert f"{name}={repr(getattr(snapshot, name))}" in detail, key
        assert f"s={repr(snapshot.active_end)}" in detail, key
        assert snapshot.solver_route in detail, key
        for turn in snapshot.turns:
            assert repr(turn.x) in detail, key
        for candidate in snapshot.candidates:
            assert f"Slot {candidate.slot_idx}" in detail, key
            assert candidate.reason in detail, key
            if candidate.physical:
                assert repr(candidate.x) in detail, key


def test_summary_separates_geometric_turns_from_packet_interval_turns() -> None:
    """A near-endpoint reversal is visible without claiming a packet split."""
    b = 1.2
    x_0 = -math.sqrt(9 - b**2)
    negative_turn_x = x_0 + 0.25 * CLOSE_LINE_THRESHOLD * 2
    q = -((negative_turn_x**2 + b**2) ** 1.5) / b**2
    inputs = ExplorerInputs(1, 3, 1 + q, 3 + q, b, x_0, 0)
    application = ResonanceExplorerWidget(
        initial_inputs=inputs,
        renderer=lambda _analysis, _selection: Figure(),
    )

    assert len(application.analysis.turns) == 2
    assert len(application.analysis.transport_turns) == 0
    assert "Active geometric turns: 2" in application.summary_pane.object
    assert "packet interval turns: 0" in application.summary_pane.object


def test_summary_explains_target_miss_against_active_velocity_range() -> None:
    """A target above the active curve has an explicit reader-facing reason."""
    inputs = replace(PRESETS[DEFAULT_PRESET].inputs, w=0.7)
    application = ResonanceExplorerWidget(
        initial_inputs=inputs,
        renderer=lambda _analysis, _selection: Figure(),
    )

    assert application.analysis.crossings == ()
    assert "Active projected-velocity range" in application.summary_pane.object
    assert "outside active projected-velocity range" in (
        application.summary_pane.object
    )


def test_table_displays_analytic_and_packet_gradients_separately() -> None:
    """Near-zero gradients retain the packet-facing exception predicate."""
    inputs = ExplorerInputs(1, 3, 1e-13, 3e-13, 0, 1, 2e-13)
    application = ResonanceExplorerWidget(
        initial_inputs=inputs,
        renderer=lambda _analysis, _selection: Figure(),
    )

    frame = application.candidate_table.value
    assert "G" in frame
    assert "packet G" in frame
    assert any(frame["packet G"].fillna(0) != 0)


def test_candidate_row_selection_is_read_only_and_highlights_raw_slot() -> None:
    """Clicking a candidate selects its raw slot without entering cell edit."""
    rendered_slots: list[int | None] = []

    def render_selection(
        _analysis: ExampleAnalysis, slot_idx: int | None
    ) -> Figure:
        rendered_slots.append(slot_idx)
        return Figure()

    application = ResonanceExplorerWidget(renderer=render_selection)
    assert application.candidate_table.disabled
    assert application.candidate_table.selectable == 1

    application.candidate_table.selection = [0]

    selected_slot = int(application.candidate_table.value.iloc[0]["slot"])
    assert application.selected_slot_index == selected_slot
    assert rendered_slots == [None, selected_slot]


def test_guided_homology_is_default_with_seven_navigable_lessons() -> None:
    """The first screen teaches the controlled comparison before the catalog."""
    application = ResonanceExplorerWidget(renderer=lambda _analysis, _selection: Figure())

    assert len(LESSONS) == 7
    assert application.lesson_select.value == "homology"
    assert application.analysis.inputs.v_in == 1
    assert application.analysis.inputs.v_out == 3
    assert application.analysis.solver_route == "zero-intercept Q=0"
    assert application.comparison_toggle.visible
    assert application.reference_panels.active == []
    assert "v(r_in)=1" in application.lesson_note.object
    assert "v(r_out)=3" in application.lesson_note.object


def test_switching_homology_and_linked_changes_only_velocity_intercept() -> (
    None
):
    """Geometry and target stay fixed while the comparison law changes."""
    application = ResonanceExplorerWidget(renderer=lambda _analysis, _selection: Figure())
    reference = application.analysis.inputs

    application.lesson_select.value = "linked"

    linked = application.analysis.inputs
    assert (linked.r_in, linked.r_out, linked.b, linked.x_0, linked.w) == (
        reference.r_in,
        reference.r_out,
        reference.b,
        reference.x_0,
        reference.w,
    )
    assert linked.v_in == -1.34375
    assert linked.v_out == 0.65625
    assert len(application.analysis.crossings) == 3
    assert "Original equation" in application.lesson_note.object
    assert "Squared solver equation" in application.lesson_note.object


def test_guided_pair_states_change_only_the_named_causal_input() -> None:
    """Launch and target cards retain the linked law and other coordinates."""
    application = ResonanceExplorerWidget(renderer=lambda _analysis, _selection: Figure())
    application.lesson_select.value = "linked"
    linked = application.analysis.inputs

    application.lesson_select.value = "launch"
    application.lesson_state_select.value = "two_crossings"
    two = application.analysis.inputs
    assert two.x_0 == -0.5
    assert replace(two, x_0=linked.x_0) == linked
    application.lesson_state_select.value = "one_crossing"
    one = application.analysis.inputs
    assert one.x_0 == 1.2
    assert replace(one, x_0=linked.x_0) == linked

    application.lesson_select.value = "target"
    application.lesson_state_select.value = "above_range"
    miss = application.analysis.inputs
    assert miss.w == 0.7
    assert replace(miss, w=linked.w) == linked


def test_every_guided_state_is_reachable_without_opening_reference_panels() -> (
    None
):
    """The seven-card route covers every paired state and current solver route."""
    application = ResonanceExplorerWidget(renderer=lambda _analysis, _selection: Figure())

    for lesson in LESSONS:
        application.lesson_select.value = lesson.key
        for state in lesson.states:
            application.lesson_state_select.value = state.key
            assert application.analysis.inputs == state.inputs
            assert application.analysis.solver_route == state.expected_route
            assert len(application.slot_chips.options) == 4
            assert application.reference_panels.active == []


def test_slot_chips_show_four_raw_decisions_and_select_slot() -> None:
    """The linked fourth slot is exposed as a wrong-sign artifact."""
    rendered_slots: list[int | None] = []

    def render_selection(_analysis: object, slot_idx: int | None) -> Figure:
        rendered_slots.append(slot_idx)
        return Figure()

    application = ResonanceExplorerWidget(renderer=render_selection)
    application.lesson_select.value = "linked"
    options = application.slot_chips.options

    assert len(options) == 4
    assert any("wrong sign" in label for label in options)
    wrong_slot = next(
        candidate.slot_idx
        for candidate in application.analysis.candidates
        if "squared-sign artifact" in candidate.flags
    )
    application.slot_chips.value = wrong_slot

    assert application.selected_slot_index == wrong_slot
    assert "wrong unsquared sign" in application.slot_fate_pane.object
    assert rendered_slots[-1] == wrong_slot
    assert "4/4 real slots" in application.pipeline_pane.object
    assert "3 original-equation matches" in application.pipeline_pane.object


def test_custom_edit_keeps_reference_panels_available() -> None:
    """Exact editing leaves the user in a computed custom state."""
    application = ResonanceExplorerWidget(renderer=lambda _analysis, _selection: Figure())

    application.controls["w"].value = 0.4

    assert application.lesson_select.value is None
    assert application.lesson_state_select.disabled
    assert application.analysis.inputs.w == 0.4
    assert application.figure_pane.object is not None
    assert len(application.candidate_table.value) == 4
    assert application.inputs_panel in application.reference_panels.objects


def test_catalog_search_filters_without_replacing_current_case() -> None:
    """Searching the full catalog is read-only until a named case is chosen."""
    application = ResonanceExplorerWidget(renderer=lambda _analysis, _selection: Figure())
    starting_inputs = application.analysis.inputs

    application.catalog_search.value = "radial inward"

    assert len(application.preset_select.options) == 2
    assert "radial inward" in application.preset_select.options.values()
    assert application.analysis.inputs == starting_inputs


def test_evidence_panel_clears_on_invalid_input_and_tracks_new_snapshot(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Independent evidence cannot survive an input edit as stale output."""
    observed_targets: list[float] = []

    def fake_evidence(analysis: object) -> dict[str, object]:
        observed_targets.append(analysis.inputs.w)
        return {
            "inputs": {"w": analysis.inputs.w},
            "agreement": {"status": "matched"},
            "packet_nearest": {"agreement_within_1e-8": True},
        }

    monkeypatch.setattr(
        import_module(ResonanceExplorerWidget.__module__),
        "build_evidence_record",
        fake_evidence,
    )
    application = ResonanceExplorerWidget(renderer=lambda _analysis, _selection: Figure())

    application._on_evidence_request(None)
    assert application.evidence_pane.object["inputs"]["w"] == 0.3

    application.controls["w"].value = None
    assert application.evidence_pane.object == {}
    assert application.evidence_button.disabled

    application.controls["w"].value = 0.4
    assert application.evidence_pane.object == {}
    assert not application.evidence_button.disabled
    application._on_evidence_request(None)
    assert observed_targets == [0.3, 0.4]
    assert application.evidence_pane.object["inputs"]["w"] == 0.4
