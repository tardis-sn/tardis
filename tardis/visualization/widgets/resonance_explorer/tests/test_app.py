"""Orchestration checks for the local explorer interface."""

# ruff: noqa: TID252 - Package-local tests use relative imports.

from __future__ import annotations

from dataclasses import replace

from bokeh.document import Document
from matplotlib.figure import Figure

from ..app import ResonanceExplorerWidget
from ..lessons import LESSONS


def test_default_case_serializes_to_a_browser_document() -> None:
    """The guided layout must contain only values Panel can serialize."""
    application = ResonanceExplorerWidget()
    document = Document()

    application.layout.server_doc(document)

    assert document.roots


def test_guided_homology_is_default_with_seven_navigable_lessons() -> None:
    """The first screen teaches the controlled comparison."""
    application = ResonanceExplorerWidget(renderer=lambda _analysis, _selection: Figure())

    assert len(LESSONS) == 7
    assert application.lesson_select.value == "homology"
    assert application.analysis.inputs.v_in == 1
    assert application.analysis.inputs.v_out == 3
    assert application.analysis.solver_route == "zero-intercept Q=0"
    assert application.comparison_toggle.visible
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


def test_every_guided_state_is_reachable() -> None:
    """The seven-card route covers every paired state and current solver route."""
    application = ResonanceExplorerWidget(renderer=lambda _analysis, _selection: Figure())

    for lesson in LESSONS:
        application.lesson_select.value = lesson.key
        for state in lesson.states:
            application.lesson_state_select.value = state.key
            assert application.analysis.inputs == state.inputs
            assert application.analysis.solver_route == state.expected_route
            assert len(application.slot_chips.options) == 4


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
