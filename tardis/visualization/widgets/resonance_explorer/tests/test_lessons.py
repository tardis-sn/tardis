"""Checks that each lesson step takes the solver route its text states."""

import pytest

from tardis.visualization.widgets.resonance_explorer.analysis import (
    analyze_case,
)
from tardis.visualization.widgets.resonance_explorer.lessons import (
    LESSONS,
    get_lesson_state,
)


@pytest.mark.parametrize(
    ("lesson_key", "state_key"),
    [(lesson.key, state.key) for lesson in LESSONS for state in lesson.states],
)
def test_lesson_state_takes_its_stated_route(
    lesson_key: str, state_key: str
) -> None:
    """Claim: the route named in the lesson text is the computed route."""
    state = get_lesson_state(lesson_key, state_key)
    assert analyze_case(state.inputs).solver_route == state.expected_route
