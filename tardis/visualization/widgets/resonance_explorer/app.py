"""Local Panel controls and linked views for single-shell resonance."""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

import panel as pn

from .analysis import CaseAnalysis, analyze_case
from .lessons import LESSONS, LESSONS_BY_KEY, get_lesson_state
from .model import ExplorerInputs
from .render import render_guided_case

if TYPE_CHECKING:
    from matplotlib.figure import Figure


def _format_number(value: float | None) -> str:
    """Show an exact zero while making an absent value explicit."""
    return "none" if value is None else f"{value:.9g}"


def _lesson_options() -> dict[str, str | None]:
    """Present the short causal sequence of guided lessons."""
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


class ResonanceExplorerWidget:
    """Own guided lesson controls and one committed analysis snapshot.

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
        pn.extension()
        self.analyzer = analyzer or analyze_case
        self.renderer = renderer
        self.guided_renderer = (
            None if renderer is not None else render_guided_case
        )
        self.analysis: CaseAnalysis | None = None
        self.selected_slot_index: int | None = None
        self._setting_lesson = False

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
        self.figure_pane = pn.pane.Matplotlib(
            None, tight=True, sizing_mode="stretch_width"
        )
        self.layout = self._build_layout()

        self.lesson_select.param.watch(self._on_lesson_change, "value")
        self.lesson_state_select.param.watch(
            self._on_lesson_state_change, "value"
        )
        self.comparison_toggle.param.watch(self._on_comparison_change, "value")
        self.slot_chips.param.watch(self._on_slot_chip_change, "value")
        if initial_inputs is None:
            self._load_lesson("homology", "initial")
        else:
            self.set_inputs(initial_inputs)

    def _build_layout(self) -> pn.layout.Row:
        """Lead with the short lesson sequence beside the linked plot."""
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
        primary = pn.Column(
            self.headline_pane,
            self.decision_pane,
            self.figure_pane,
            "### Four production-ordered root slots",
            self.slot_chips,
            self.slot_fate_pane,
            self.pipeline_pane,
            sizing_mode="stretch_width",
        )
        return pn.Row(
            navigation,
            primary,
            sizing_mode="stretch_width",
        )

    def _mark_guided_custom(self) -> None:
        """Identify a non-lesson case without discarding its computed snapshot."""
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
            self.lesson_note.object = "Choose a guided lesson."
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

    def _commit_inputs(self, inputs: ExplorerInputs) -> None:
        """Compute a complete snapshot before changing any scientific view."""
        snapshot = self.analyzer(inputs)
        figure = self._render_snapshot(snapshot, None)
        old_figure = self.figure_pane.object
        self.analysis = snapshot
        self.selected_slot_index = None
        self.figure_pane.object = figure
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
        if old_figure is not None:
            old_figure.clear()

    def set_inputs(self, inputs: ExplorerInputs) -> None:
        """Commit exact shell, ray, target, and stop values in one update.

        Parameters
        ----------
        inputs : ExplorerInputs
            Canonical values to analyze.
        """
        if not self._setting_lesson:
            self._mark_guided_custom()
        self._commit_inputs(inputs)

    def display(self) -> pn.layout.Row:
        """Return the live Panel layout for a notebook or local server.

        Returns
        -------
        panel.layout.Row
            Lesson controls, linked plot, and decisions.
        """
        return self.layout


def build_app() -> pn.layout.Row:
    """Create an independent browser-session explorer layout.

    Returns
    -------
    panel.layout.Row
        Locally served lesson controls, linked figure, and decisions.
    """
    return ResonanceExplorerWidget().display()
