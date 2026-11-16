"""Local Panel controls and linked views for single-shell resonance."""

from __future__ import annotations

import io
import math
from collections.abc import Callable
from dataclasses import fields
from typing import TYPE_CHECKING

import pandas as pd
import panel as pn

from .analysis import CaseAnalysis, analyze_case
from .diagnostics import (
    GENERIC_SOLVER_EXAMPLES,
    diagnose_generic_quartic,
    diagnose_solver,
)
from .evidence import build_evidence_record
from .export import export_figure_bytes, export_json_bytes, json_safe
from .fingerprint import input_fingerprint
from .geometry import calculate_projected_gradient
from .lessons import LESSONS, LESSONS_BY_KEY, get_lesson_state
from .model import ExplorerInputs
from .presets import DEFAULT_PRESET, PRESETS
from .render import render_case, render_guided_case

if TYPE_CHECKING:
    from matplotlib.figure import Figure


INPUT_LABELS = {
    "r_in": "Inner radius r_in",
    "r_out": "Outer radius r_out",
    "v_in": "Radial velocity at r_in",
    "v_out": "Radial velocity at r_out",
    "b": "Impact parameter b",
    "x_0": "Signed launch coordinate x_0",
    "w": "Line target w",
}
TABLE_COLUMNS = (
    "slot",
    "X",
    "x",
    "s",
    "r",
    "G",
    "packet G",
    "residual",
    "unsquared tolerance",
    "production accepted",
    "production event eligible",
    "unsquared verified",
    "verified event eligible",
    "multiplicity / contact",
    "flags",
    "reason",
)


def _format_number(value: float | None) -> str:
    """Show an exact zero while making an absent value explicit."""
    return "none" if value is None else f"{value:.9g}"


def _preset_options() -> dict[str, str | None]:
    """Order examples by nominal ray, projected shape, and line outcome."""
    options = {}
    for key, preset in PRESETS.items():
        outcome = (
            "contact or limit"
            if preset.expected_crossing_count is None
            else f"{preset.expected_crossing_count} crossings"
        )
        label = (
            f"{preset.expected_ray_family} → {preset.expected_shape_family} "
            f"→ {outcome} · {preset.name}"
        )
        options[label] = key
    options["Custom / modified inputs"] = None
    return options


def _lesson_options() -> dict[str, str | None]:
    """Present the short causal sequence before the full reference catalog."""
    return {
        **{lesson.title: lesson.key for lesson in LESSONS},
        "Custom / modified inputs": None,
    }


def _short_slot_fate(candidate: object) -> str:
    """Summarize a raw slot without conflating solver and packet decisions."""
    flags = candidate.flags
    if candidate.x is None:
        return "non-real"
    if "squared-sign artifact" in flags:
        return "wrong unsquared sign"
    if not candidate.unsquared_verified:
        return "equation-unverified"
    if any(
        flag in flags
        for flag in (
            "at or behind start",
            "after event stop",
            "outside shell",
            "past strict caller end",
            "boundary tie",
        )
    ):
        return "outside active flight"
    return "eligible" if candidate.verified_event_eligible else candidate.reason


def _slot_chip_label(candidate: object) -> str:
    """Keep all four choices readable at a typical browser width."""
    short = {
        "wrong unsquared sign": "wrong sign",
        "outside active flight": "outside",
        "equation-unverified": "unverified",
    }.get(_short_slot_fate(candidate), _short_slot_fate(candidate))
    return f"{candidate.slot_idx} · {short}"


def _guided_headline(analysis: CaseAnalysis) -> str:
    """State the causal chain in one reader-facing sentence."""
    return (
        f"**{analysis.ray_family} → {analysis.shape_family} → "
        f"{analysis.result_family}.**"
    )


def _guided_decision(analysis: CaseAnalysis) -> str:
    """Report route, flight event, and potential resonance independently."""
    nearest = analysis.nearest_resonance
    production = getattr(analysis, "production_nearest_resonance", nearest)
    verified_text = (
        "none"
        if nearest is None
        else f"slot {nearest.slot_idx} at s={_format_number(nearest.s)}"
    )
    production_text = (
        "none"
        if production is None
        else f"slot {production.slot_idx} at s={_format_number(production.s)}"
    )
    return (
        f"**Production route:** {analysis.solver_route} · "
        f"**First boundary:** {analysis.boundary_kind} at "
        f"s={_format_number(analysis.boundary_distance)} · "
        f"**Active end:** {analysis.stop_reason} at "
        f"s={_format_number(analysis.active_end)}  \n"
        f"**Potential resonance:** nearest equation-verified eligible "
        f"{verified_text}; production-selected {production_text}. "
        "Absorption requires opacity and sampled optical depth."
    )


def _guided_pipeline(analysis: CaseAnalysis) -> str:
    """Count each filtering stage using the committed analysis snapshot."""
    real_count = sum(
        candidate.x is not None for candidate in analysis.candidates
    )
    verified_count = sum(
        candidate.unsquared_verified for candidate in analysis.candidates
    )
    eligible_count = sum(
        candidate.verified_event_eligible for candidate in analysis.candidates
    )
    next_event = (
        f"potential resonance at s={_format_number(analysis.nearest_resonance.s)}"
        if analysis.nearest_resonance is not None
        else f"{analysis.stop_reason} at s={_format_number(analysis.active_end)}"
    )
    return (
        f"**Decision path:** {analysis.solver_route} → "
        f"{real_count}/4 real slots → {verified_count} original-equation "
        f"matches → {eligible_count} verified forward candidates → {next_event}."
    )


def _candidate_frame(analysis: CaseAnalysis) -> pd.DataFrame:
    """Keep every production slot in its original order."""
    finite_coordinates = [
        candidate.x
        for candidate in analysis.candidates
        if candidate.x is not None and math.isfinite(candidate.x)
    ]
    contacts = getattr(analysis, "analytic_contacts", ())

    def contact_status(candidate: object) -> str:
        """Describe exact slot duplication separately from geometric contact."""
        if candidate.x is None or not math.isfinite(candidate.x):
            return "unavailable"
        repetitions = finite_coordinates.count(candidate.x)
        labels = [
            f"{repetitions} identical raw slots"
            if repetitions > 1
            else "one raw slot"
        ]
        for contact in contacts:
            if candidate.x == contact.x:
                labels.append(f"analytic {contact.kind} at this x")
            elif math.isclose(candidate.x, contact.x, rel_tol=0, abs_tol=1e-9):
                labels.append(
                    f"near analytic {contact.kind}; multiplicity unconfirmed"
                )
        return "; ".join(labels)

    rows = [
        {
            "slot": candidate.slot_idx,
            "X": candidate.scaled_root,
            "x": candidate.x,
            "s": candidate.s,
            "r": candidate.r,
            "G": candidate.gradient,
            "packet G": candidate.packet_gradient,
            "residual": candidate.residual,
            "unsquared tolerance": candidate.unsquared_tolerance,
            "production accepted": candidate.accepted,
            "production event eligible": candidate.event_eligible,
            "unsquared verified": candidate.unsquared_verified,
            "verified event eligible": candidate.verified_event_eligible,
            "multiplicity / contact": contact_status(candidate),
            "flags": ", ".join(candidate.flags),
            "reason": candidate.reason,
        }
        for candidate in analysis.candidates
    ]
    return pd.DataFrame(rows, columns=TABLE_COLUMNS)


def _case_summary(analysis: CaseAnalysis) -> str:
    """Explain the case and distinguish resonance from interaction."""
    inputs = analysis.inputs
    nearest = analysis.nearest_resonance
    production_nearest = getattr(
        analysis, "production_nearest_resonance", nearest
    )
    unverified_accepted = sum(
        candidate.accepted and not candidate.unsquared_verified
        for candidate in analysis.candidates
    )
    nearest_text = (
        "none"
        if nearest is None
        else f"slot {nearest.slot_idx} at x={_format_number(nearest.x)}, "
        f"s={_format_number(nearest.s)}"
    )
    production_nearest_text = (
        "none"
        if production_nearest is None
        else f"slot {production_nearest.slot_idx} at "
        f"x={_format_number(production_nearest.x)}, "
        f"s={_format_number(production_nearest.s)} "
        f"(unsquared {'verified' if production_nearest.unsquared_verified else 'unverified'})"
    )
    if nearest is None:
        next_event = f"{analysis.stop_reason} first at s={_format_number(analysis.active_end)}"
    elif nearest.s == analysis.active_end:
        next_event = (
            f"potential line at s={_format_number(nearest.s)} ties the assumed stop; "
            "the strict scatter comparison checks the line candidate first"
        )
    else:
        next_event = (
            f"potential line at s={_format_number(nearest.s)} before "
            f"{analysis.stop_reason}"
        )
    lines = [
        "### Current case",
        f"**{analysis.ray_family} → {analysis.shape_family} → {analysis.result_family}**",
        f"First boundary: **{analysis.boundary_kind}** at "
        f"s={_format_number(analysis.boundary_distance)}. "
        f"Active flight ends at {analysis.stop_reason}, s={_format_number(analysis.active_end)}.",
        f"Radial law: **{analysis.radial_law}**, "
        f"v(r)={_format_number(analysis.m)} r + ({_format_number(analysis.q)}). "
        f"Solver route: **{analysis.solver_route}**.",
        f"Launch: b={_format_number(inputs.b)}, x_0={_format_number(inputs.x_0)}, "
        f"r_0={_format_number(analysis.r_0)}, μ_0={_format_number(analysis.mu_0)}. "
        f"Shell: [{_format_number(inputs.r_in)}, {_format_number(inputs.r_out)}]. "
        f"Edge velocities: [{_format_number(inputs.v_in)}, {_format_number(inputs.v_out)}].",
        f"Target w={_format_number(inputs.w)} = c(1 − ν_line/ν_packet). "
        f"Active geometric turns: {len(analysis.turns)}; packet interval turns: "
        f"{len(getattr(analysis, 'transport_turns', analysis.turns))}; "
        "monotonic intervals: "
        f"{len(getattr(analysis, 'monotonic_intervals', ()))}; distinct verified "
        f"forward crossings: {len(analysis.crossings)}.",
        f"Production accepted but unsquared unverified: {unverified_accepted}. "
        "These raw slots are excluded from verified crossings and verified "
        "event selection.",
        f"Production-selected line candidate: **{production_nearest_text}**. "
        "This reports the current transport decision separately from "
        "unsquared-equation verification.",
        f"Nearest verified eligible resonance: **{nearest_text}**. A resonance is a "
        "potential line encounter; actual absorption needs opacity and sampled optical depth.",
        f"Next geometric opportunity: **{next_event}**.",
    ]
    velocity_range = getattr(analysis, "active_velocity_range", None)
    if velocity_range is not None:
        lines.append(
            "Active projected-velocity range: "
            f"[{_format_number(velocity_range[0])}, "
            f"{_format_number(velocity_range[1])}]."
        )
    no_crossing_reason = getattr(analysis, "no_crossing_reason", None)
    if no_crossing_reason:
        lines.append(f"Why no verified crossing: {no_crossing_reason}")
    if getattr(analysis, "notes", ()):
        lines.extend(f"- {note}" for note in analysis.notes)
    for interval_idx, (x_start, x_end) in enumerate(
        getattr(analysis, "monotonic_intervals", ()), start=1
    ):
        midpoint = (x_start + x_end) / 2
        gradient = calculate_projected_gradient(midpoint, inputs)
        direction = (
            "decreases"
            if gradient > 0
            else "increases"
            if gradient < 0
            else "flat"
        )
        lines.append(
            f"Interval {interval_idx}: x from {_format_number(x_start)} to "
            f"{_format_number(x_end)}; comoving frequency {direction} "
            "as s increases, since dν_comoving/ds = −(ν_packet/c)G."
        )
    return "\n\n".join(lines)


def _case_detail(analysis: CaseAnalysis, preset_key: str | None) -> str:
    """Describe exact named or custom inputs using the committed calculation."""
    inputs = analysis.inputs
    if preset_key in PRESETS:
        preset = PRESETS[preset_key]
        heading = f"**{preset.name}.** {preset.description}"
        crosswalk = f"Manuscript crosswalk: {preset.source_crosswalk}."
    else:
        heading = "**Custom / modified inputs.**"
        crosswalk = "The current values are not a named preset."
    turns = (
        ", ".join(
            f"x={repr(turn.x)} (s={repr(turn.s)}, {turn.kind})"
            for turn in analysis.turns
        )
        or "none"
    )
    physical_roots = (
        ", ".join(
            f"slot {candidate.slot_idx}: x={repr(candidate.x)} (s={repr(candidate.s)})"
            for candidate in analysis.candidates
            if candidate.physical
        )
        or "none"
    )
    lines = [
        heading,
        crosswalk,
        "**Exact inputs and derived radial law**",
        f"r_in={repr(inputs.r_in)}; r_out={repr(inputs.r_out)}; "
        f"v_in={repr(inputs.v_in)}; v_out={repr(inputs.v_out)}.",
        f"b={repr(inputs.b)}; x_0={repr(inputs.x_0)}; "
        f"w={repr(inputs.w)}; s_stop={repr(inputs.s_stop)}.",
        f"m={repr(analysis.m)}; q={repr(analysis.q)}.",
        f"**Active end:** s={repr(analysis.active_end)} "
        f"({analysis.stop_reason}).",
        f"**Active geometric turns ({len(analysis.turns)}):** {turns}.",
        f"**Packet interval turns:** "
        f"{len(getattr(analysis, 'transport_turns', analysis.turns))}.",
        f"**Unsquared verified physical root slots:** {physical_roots}.",
        f"**Solver route:** {analysis.solver_route}.",
        "**Candidate dispositions in raw slot order**",
    ]
    for candidate in analysis.candidates:
        flags = ", ".join(candidate.flags) or "none"
        lines.append(
            f"- Slot {candidate.slot_idx}: x={repr(candidate.x)}; "
            f"production {'accepted' if candidate.accepted else 'rejected'}; "
            f"production event {'eligible' if candidate.event_eligible else 'ineligible'}; "
            f"unsquared {'verified' if candidate.unsquared_verified else 'unverified'}; "
            f"verified event {'eligible' if candidate.verified_event_eligible else 'ineligible'}; "
            f"reason: {candidate.reason}; flags: {flags}."
        )
    if not analysis.candidates:
        lines.append("No raw quartic slots on this solver route.")
    return "\n\n".join(lines)


class ResonanceExplorerWidget:
    """Own exact controls and one committed analysis snapshot.

    Parameters
    ----------
    initial_inputs : ExplorerInputs or None, optional
        Starting values; homologous expansion opens the guided sequence.
    analyzer : callable or None, optional
        Scientific analysis boundary, injectable for orchestration tests.
    renderer : callable or None, optional
        Headless figure boundary, injectable for orchestration tests.
    """

    def __init__(
        self,
        initial_inputs: ExplorerInputs | None = None,
        analyzer: Callable[[ExplorerInputs], CaseAnalysis] | None = None,
        renderer: Callable[[CaseAnalysis, int | None], Figure] | None = None,
    ) -> None:
        pn.extension("tabulator")
        self.analyzer = analyzer or analyze_case
        self.renderer = renderer or render_case
        self.guided_renderer = (
            None if renderer is not None else render_guided_case
        )
        self.analysis: CaseAnalysis | None = None
        self.selected_slot_index: int | None = None
        self._setting_controls = False
        self._setting_lesson = False

        initial = (
            initial_inputs or get_lesson_state("homology", "initial").inputs
        )
        self.lesson_select = pn.widgets.Select(
            label="Guided lessons",
            options=_lesson_options(),
            value="homology" if initial_inputs is None else None,
            sizing_mode="stretch_width",
        )
        self.lesson_state_select = pn.widgets.RadioButtonGroup(
            label="Lesson state",
            options={},
            value=None,
            sizing_mode="stretch_width",
        )
        self.comparison_toggle = pn.widgets.Checkbox(
            label="Show homologous / nonhomologous comparison",
            value=True,
        )
        self.lesson_question = pn.pane.Markdown("", sizing_mode="stretch_width")
        self.lesson_note = pn.pane.Markdown("", sizing_mode="stretch_width")
        self.headline_pane = pn.pane.Markdown("", sizing_mode="stretch_width")
        self.decision_pane = pn.pane.Markdown("", sizing_mode="stretch_width")
        self.pipeline_pane = pn.pane.Markdown("", sizing_mode="stretch_width")
        self.slot_chips = pn.widgets.RadioButtonGroup(
            label="Production root slots",
            options={},
            value=None,
            sizing_mode="stretch_width",
        )
        self.slot_fate_pane = pn.pane.Markdown("", sizing_mode="stretch_width")
        self.controls = {
            name: pn.widgets.FloatInput(
                label=label, value=getattr(initial, name)
            )
            for name, label in INPUT_LABELS.items()
        }
        self.stop_enabled = pn.widgets.Checkbox(
            label="Assume an earlier electron scatter",
            value=initial.s_stop is not None,
        )
        self.stop_distance = pn.widgets.FloatInput(
            label="Earlier stop distance s_stop",
            value=initial.s_stop if initial.s_stop is not None else 0.5,
            disabled=initial.s_stop is None,
        )
        self.preset_select = pn.widgets.Select(
            label="Explore a named case: ray path → shape → line result",
            options=_preset_options(),
            value=DEFAULT_PRESET,
            sizing_mode="stretch_width",
        )
        self.catalog_search = pn.widgets.TextInput(
            label="Search the 31-case catalog",
            placeholder="ray, turn, route, or result",
        )
        self.preset_detail = pn.pane.Markdown("", sizing_mode="stretch_width")
        self.reset_button = pn.widgets.Button(
            label="Reset to linked three crossings"
        )
        self.error_pane = pn.pane.Alert("", alert_type="danger", visible=False)
        self.summary_pane = pn.pane.Markdown("", sizing_mode="stretch_width")
        self.figure_pane = pn.pane.Matplotlib(
            None, tight=True, sizing_mode="stretch_width"
        )
        self.candidate_table = pn.widgets.Tabulator(
            pd.DataFrame(columns=TABLE_COLUMNS),
            label="Production candidate slots (select a row to highlight)",
            selectable=1,
            disabled=True,
            show_index=False,
            height=220,
            sizing_mode="stretch_width",
        )
        self.diagnostic_pane = pn.pane.JSON(
            {}, depth=2, sizing_mode="stretch_width"
        )
        self.evidence_button = pn.widgets.Button(
            label="Build current-case evidence"
        )
        self.evidence_status = pn.pane.Markdown(
            "Compute independent unsquared roots and a mapped packet-distance "
            "comparison for this exact case. Source revision is in the JSON export."
        )
        self.evidence_pane = pn.pane.JSON(
            {}, depth=2, sizing_mode="stretch_width"
        )
        self.generic_select = pn.widgets.Select(
            label="Solver-only depressed quartic",
            options={
                "Choose an unphysical solver example": None,
                **{key: key for key in GENERIC_SOLVER_EXAMPLES},
            },
            value=None,
        )
        self.generic_pane = pn.pane.JSON(
            {}, depth=2, sizing_mode="stretch_width"
        )
        self.pdf_download = pn.widgets.FileDownload(
            label="Download PDF",
            callback=lambda: self._download_figure("pdf"),
            filename="resonance.pdf",
            color="primary",
        )
        self.png_download = pn.widgets.FileDownload(
            label="Download PNG",
            callback=lambda: self._download_figure("png"),
            filename="resonance.png",
        )
        self.json_download = pn.widgets.FileDownload(
            label="Download JSON",
            callback=self._download_json,
            filename="resonance.json",
        )
        self.layout = self._build_layout()

        self.lesson_select.param.watch(self._on_lesson_change, "value")
        self.lesson_state_select.param.watch(
            self._on_lesson_state_change, "value"
        )
        self.comparison_toggle.param.watch(self._on_comparison_change, "value")
        self.slot_chips.param.watch(self._on_slot_chip_change, "value")
        self.catalog_search.param.watch(self._on_catalog_search, "value")
        for control in self.controls.values():
            control.param.watch(self._on_control_change, "value")
        self.stop_enabled.param.watch(self._on_stop_change, "value")
        self.stop_distance.param.watch(self._on_control_change, "value")
        self.preset_select.param.watch(self._on_preset_change, "value")
        self.generic_select.param.watch(self._on_generic_change, "value")
        self.candidate_table.param.watch(self._on_selection_change, "selection")
        self.reset_button.on_click(self._on_reset)
        self.evidence_button.on_click(self._on_evidence_request)
        if initial_inputs is None:
            self._load_lesson("homology", "initial")
        else:
            self.set_inputs(initial)

    def _build_layout(self) -> pn.layout.Row:
        """Lead with the short lesson sequence and reveal research controls."""
        navigation = pn.Column(
            "## Single-shell resonance explorer",
            "Illustrative length and velocity units. Forward packet distance "
            "increases with signed ray coordinate x.",
            self.lesson_select,
            self.lesson_state_select,
            self.comparison_toggle,
            self.lesson_question,
            self.lesson_note,
            width=320,
        )
        self.inputs_panel = pn.Column(
            *self.controls.values(),
            self.stop_enabled,
            self.stop_distance,
        )
        self.catalog_panel = pn.Column(
            self.catalog_search,
            self.preset_select,
            self.preset_detail,
            self.reset_button,
        )
        self.diagnostics_panel = pn.Column(
            self.summary_pane,
            self.candidate_table,
            self.diagnostic_pane,
        )
        self.solver_examples_panel = pn.Column(
            self.generic_select, self.generic_pane
        )
        self.evidence_panel = pn.Column(
            self.evidence_button, self.evidence_status, self.evidence_pane
        )
        self.exports_panel = pn.Row(
            self.pdf_download, self.png_download, self.json_download
        )
        self.reference_panels = pn.Accordion(
            ("Exact input editing", self.inputs_panel),
            ("Searchable 31-case reference catalog", self.catalog_panel),
            (
                "Full candidate table and source-reconstructed diagnostics",
                self.diagnostics_panel,
            ),
            ("Independent guided evidence", self.evidence_panel),
            ("Exports", self.exports_panel),
            (
                "Solver-only examples (not shell trajectories)",
                self.solver_examples_panel,
            ),
            active=[],
            sizing_mode="stretch_width",
        )
        primary = pn.Column(
            self.headline_pane,
            self.decision_pane,
            self.error_pane,
            self.figure_pane,
            "### Four production-ordered root slots",
            self.slot_chips,
            self.slot_fate_pane,
            self.pipeline_pane,
            self.reference_panels,
            sizing_mode="stretch_width",
        )
        return pn.Row(
            navigation,
            primary,
            sizing_mode="stretch_width",
        )

    def _mark_guided_custom(self) -> None:
        """Identify an edited case without discarding its computed snapshot."""
        if self._setting_lesson:
            return
        self._setting_lesson = True
        try:
            self.lesson_select.value = None
            self.lesson_state_select.options = {}
            self.lesson_state_select.value = None
            self.lesson_state_select.disabled = True
            self.comparison_toggle.visible = False
            self.lesson_question.object = "**Custom / modified inputs.**"
            self.lesson_note.object = (
                "Inspect the exact inputs and reference catalog below, "
                "or choose a guided lesson."
            )
        finally:
            self._setting_lesson = False

    def _load_lesson(self, lesson_key: str, state_key: str) -> None:
        """Apply one exact lesson state as a single committed input change."""
        lesson = LESSONS_BY_KEY[lesson_key]
        state = get_lesson_state(lesson_key, state_key)
        self._setting_lesson = True
        try:
            self.lesson_select.value = lesson_key
            self.lesson_state_select.options = {
                item.label: item.key for item in lesson.states
            }
            self.lesson_state_select.value = state_key
            self.lesson_state_select.disabled = len(lesson.states) == 1
            self.comparison_toggle.visible = lesson_key in (
                "homology",
                "linked",
            )
            self.lesson_question.object = f"### {lesson.question}"
            self.lesson_note.object = (
                f"**Change:** {state.highlight_change}  \n"
                f"**Expected:** {state.expected_result}  \n"
                f"**Route:** {state.expected_route}  \n"
                f"{state.explanation}  \n"
                f"**Shell-edge velocities:** "
                f"v(r_in)={_format_number(state.inputs.v_in)}, "
                f"v(r_out)={_format_number(state.inputs.v_out)}."
            )
            if lesson_key == "linked":
                self.lesson_note.object += (
                    "\n\n**Original equation:** "
                    "m·x + q·x/√(x²+b²) = w  \n"
                    "**Squared solver equation:** "
                    "(m·x − w)²(x²+b²) − q²x² = 0. "
                    "Squaring loses the sign of the original equation."
                )
            self.set_inputs(state.inputs)
        finally:
            self._setting_lesson = False

    def _reference_snapshot(
        self, snapshot: CaseAnalysis
    ) -> CaseAnalysis | None:
        """Keep the controlled comparison available for fixed guided scales."""
        lesson_key = self.lesson_select.value
        if lesson_key not in ("homology", "linked"):
            return None
        paired_key = "linked" if lesson_key == "homology" else "homology"
        paired_inputs = get_lesson_state(paired_key, "initial").inputs
        if snapshot.inputs == paired_inputs:
            return None
        return self.analyzer(paired_inputs)

    def _render_snapshot(
        self, snapshot: CaseAnalysis, slot_idx: int | None
    ) -> Figure:
        """Render the guided plot while retaining injected test renderers."""
        if self.guided_renderer is None:
            return self.renderer(snapshot, slot_idx)
        reference = self._reference_snapshot(snapshot)
        return self.guided_renderer(
            snapshot,
            reference,
            slot_idx,
            show_reference=self.comparison_toggle.value,
        )

    def _on_lesson_change(self, event: object) -> None:
        """Advance to the first state of a selected lesson."""
        if self._setting_lesson or event.new not in LESSONS_BY_KEY:
            return
        lesson = LESSONS_BY_KEY[event.new]
        self._load_lesson(lesson.key, lesson.states[0].key)

    def _on_lesson_state_change(self, event: object) -> None:
        """Switch paired launch, target, or branch states exactly."""
        if self._setting_lesson or event.new is None:
            return
        lesson_key = self.lesson_select.value
        if lesson_key in LESSONS_BY_KEY:
            self._load_lesson(lesson_key, event.new)

    def _on_comparison_change(self, _event: object) -> None:
        """Redraw the same committed state with or without its reference."""
        if self.analysis is None or self.guided_renderer is None:
            return
        figure = self._render_snapshot(self.analysis, self.selected_slot_index)
        previous = self.figure_pane.object
        self.figure_pane.object = figure
        if previous is not None:
            previous.clear()

    def _on_catalog_search(self, event: object) -> None:
        """Filter named cases without changing the active scientific case."""
        query = event.new.strip().casefold()
        options = {
            label: key
            for label, key in _preset_options().items()
            if key is None
            or not query
            or any(
                query in value.casefold()
                for value in (
                    label,
                    PRESETS[key].description,
                    PRESETS[key].source_crosswalk,
                )
            )
        }
        selected = self.preset_select.value
        if selected is not None and selected not in options.values():
            self.preset_select.value = None
        self.preset_select.options = options

    def _on_slot_chip_change(self, event: object) -> None:
        """Describe and highlight the selected production-order raw slot."""
        if self.analysis is None or event.new is None:
            return
        candidate = next(
            (
                item
                for item in self.analysis.candidates
                if item.slot_idx == event.new
            ),
            None,
        )
        if candidate is None:
            return
        fate = _short_slot_fate(candidate)
        self.slot_fate_pane.object = (
            f"**Slot {candidate.slot_idx}: {fate}.** "
            f"Production {'accepted' if candidate.accepted else 'rejected'}; "
            f"original equation {'verified' if candidate.unsquared_verified else 'unverified'}; "
            f"{candidate.reason}."
        )
        figure = self._render_snapshot(self.analysis, candidate.slot_idx)
        previous = self.figure_pane.object
        self.figure_pane.object = figure
        self.selected_slot_index = candidate.slot_idx
        if previous is not None:
            previous.clear()

    def _on_evidence_request(self, _event: object) -> None:
        """Build independent evidence only for the current valid snapshot."""
        snapshot = self.analysis
        if snapshot is None or not isinstance(snapshot, CaseAnalysis):
            return
        record = build_evidence_record(snapshot)
        self.evidence_pane.object = json_safe(record)
        packet = record["packet_nearest"]
        agreement = record["agreement"]
        if agreement["status"] != "matched":
            status = "**Equation-root mismatch:** independent root missing from verified slots."
        elif packet.get("agreement_within_1e-8") is False:
            status = "**Packet-distance mismatch:** inspect the two reported distances."
        else:
            status = "Independent equation roots and production slots are reported below."
        self.evidence_status.object = status

    def _read_inputs(self) -> ExplorerInputs:
        """Read exact numeric fields and optional earlier stop."""
        values = {name: widget.value for name, widget in self.controls.items()}
        if self.stop_enabled.value:
            values["s_stop"] = self.stop_distance.value
        for name, value in values.items():
            try:
                finite = math.isfinite(value)
            except TypeError, ValueError:
                finite = False
            if not finite:
                label = INPUT_LABELS.get(name, "Earlier stop distance s_stop")
                raise ValueError(f"{label} must be a finite number.")
        values.setdefault("s_stop", None)
        return ExplorerInputs(**values)

    def _clear_numerical_views(self, message: str) -> None:
        """Remove prior scientific output after a rejected input."""
        self.analysis = None
        self.selected_slot_index = None
        previous = self.figure_pane.object
        self.figure_pane.object = None
        if previous is not None:
            previous.clear()
        self.summary_pane.object = "No current analysis: correct the input."
        self.headline_pane.object = "No current analysis: correct the input."
        self.decision_pane.object = ""
        self.pipeline_pane.object = ""
        self.slot_chips.options = {}
        self.slot_chips.value = None
        self.slot_fate_pane.object = ""
        self.evidence_pane.object = {}
        self.evidence_status.object = (
            "Correct the input before building evidence."
        )
        self.evidence_button.disabled = True
        self.candidate_table.value = pd.DataFrame(columns=TABLE_COLUMNS)
        self.diagnostic_pane.object = {}
        self.preset_detail.object = (
            "**Custom / modified inputs.** Correct the input to see "
            "computed case details."
        )
        self.error_pane.object = message
        self.error_pane.visible = True
        for download in (
            self.pdf_download,
            self.png_download,
            self.json_download,
        ):
            download.disabled = True

    def _commit_inputs(self, inputs: ExplorerInputs) -> None:
        """Compute a complete snapshot before changing any scientific view."""
        try:
            snapshot = self.analyzer(inputs)
        except ValueError as error:
            self._clear_numerical_views(str(error))
            return
        figure = self._render_snapshot(snapshot, None)
        frame = _candidate_frame(snapshot)
        summary = _case_summary(snapshot)
        detail = _case_detail(snapshot, self.preset_select.value)
        diagnostic = (
            diagnose_solver(snapshot)
            if isinstance(snapshot, CaseAnalysis)
            else {"route": snapshot.solver_route, "provenance": "test snapshot"}
        )
        old_figure = self.figure_pane.object
        self.analysis = snapshot
        self.selected_slot_index = None
        self.figure_pane.object = figure
        self.summary_pane.object = summary
        self.headline_pane.object = _guided_headline(snapshot)
        self.decision_pane.object = _guided_decision(snapshot)
        self.pipeline_pane.object = _guided_pipeline(snapshot)
        self.slot_chips.options = {
            _slot_chip_label(candidate): candidate.slot_idx
            for candidate in snapshot.candidates
        }
        self.slot_chips.value = None
        self.slot_fate_pane.object = (
            "Select a raw slot to see its disposition. Four algebraic slots "
            "do not imply four packet encounters."
        )
        self.evidence_pane.object = {}
        self.evidence_status.object = (
            "Compute independent unsquared roots and a mapped packet-distance "
            "comparison for this exact case. Source revision is in the JSON export."
        )
        self.evidence_button.disabled = not isinstance(snapshot, CaseAnalysis)
        self.preset_detail.object = detail
        self.candidate_table.value = frame
        self.candidate_table.selection = []
        self.diagnostic_pane.object = json_safe(diagnostic)
        self.error_pane.visible = False
        fingerprint = input_fingerprint(snapshot.inputs)
        for extension, download in (
            ("pdf", self.pdf_download),
            ("png", self.png_download),
            ("json", self.json_download),
        ):
            download.filename = f"resonance-{fingerprint}.{extension}"
            download.disabled = False
        if old_figure is not None:
            old_figure.clear()

    def set_inputs(self, inputs: ExplorerInputs) -> None:
        """Commit exact shell, ray, target, and stop values in one update.

        Parameters
        ----------
        inputs : ExplorerInputs
            Canonical values to display and analyze.
        """
        selected_key = self.preset_select.value
        if (
            selected_key not in PRESETS
            or inputs != PRESETS[selected_key].inputs
        ):
            self.preset_select.value = None
        if not self._setting_lesson:
            self._mark_guided_custom()
        self._setting_controls = True
        try:
            for field in fields(inputs):
                if field.name in self.controls:
                    self.controls[field.name].value = getattr(
                        inputs, field.name
                    )
            self.stop_enabled.value = inputs.s_stop is not None
            self.stop_distance.value = (
                inputs.s_stop if inputs.s_stop is not None else 0.5
            )
            self.stop_distance.disabled = inputs.s_stop is None
        finally:
            self._setting_controls = False
        self._commit_inputs(inputs)

    def _on_control_change(self, _event: object) -> None:
        """Commit one exact numeric field edit."""
        if self._setting_controls:
            return
        self.preset_select.value = None
        self._mark_guided_custom()
        try:
            inputs = self._read_inputs()
        except ValueError as error:
            self._clear_numerical_views(str(error))
            return
        self._commit_inputs(inputs)

    def _on_stop_change(self, event: object) -> None:
        """Toggle the named earlier event and its input field."""
        self.stop_distance.disabled = not self.stop_enabled.value
        self._on_control_change(event)

    def _on_preset_change(self, event: object) -> None:
        """Load every canonical field from one named example."""
        key = event.new
        if key in PRESETS:
            self._mark_guided_custom()
            self.set_inputs(PRESETS[key].inputs)

    def _on_generic_change(self, event: object) -> None:
        """Show a separately labeled nonphysical solver-only diagnostic."""
        key = event.new
        diagnostic = (
            diagnose_generic_quartic(GENERIC_SOLVER_EXAMPLES[key])
            if key in GENERIC_SOLVER_EXAMPLES
            else {}
        )
        self.generic_pane.object = json_safe(diagnostic)

    def _on_reset(self, _event: object) -> None:
        """Recover the linked three-crossing example exactly."""
        if self.preset_select.value == DEFAULT_PRESET:
            self.set_inputs(PRESETS[DEFAULT_PRESET].inputs)
        else:
            self.preset_select.value = DEFAULT_PRESET

    def _on_selection_change(self, event: object) -> None:
        """Keep the detailed table selection in sync with the slot chips."""
        if self.analysis is None or not event.new:
            return
        row_idx = event.new[0]
        frame = self.candidate_table.value
        if row_idx >= len(frame):
            return
        slot_idx = int(frame.iloc[row_idx]["slot"])
        self.slot_chips.value = slot_idx

    def _download_figure(self, format_name: str) -> io.BytesIO:
        """Capture one committed snapshot for a vector or raster download."""
        snapshot = self.analysis
        if snapshot is None:
            raise ValueError("No valid analysis is available for download")
        return io.BytesIO(
            export_figure_bytes(snapshot, format_name, self.selected_slot_index)
        )

    def _download_json(self) -> io.BytesIO:
        """Capture one committed snapshot for a strict JSON download."""
        snapshot = self.analysis
        if snapshot is None:
            raise ValueError("No valid analysis is available for download")
        return io.BytesIO(export_json_bytes(snapshot))

    def display(self) -> pn.layout.Row:
        """Return the live Panel layout for a notebook or local server.

        Returns
        -------
        panel.layout.Row
            Interactive controls, linked plots, and case details.
        """
        return self.layout


def build_app() -> pn.layout.Row:
    """Create an independent browser-session explorer layout.

    Returns
    -------
    panel.layout.Row
        Locally served controls, linked figure, decisions, and exports.
    """
    return ResonanceExplorerWidget().display()
