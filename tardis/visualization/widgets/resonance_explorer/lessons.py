"""Short, exact-input lessons for the single-shell resonance explorer."""

from dataclasses import dataclass
from math import sqrt

from .model import ExplorerInputs
from .presets import PRESETS


@dataclass(frozen=True)
class LessonState:
    """One input state and its expected observable solver behavior."""

    key: str
    label: str
    inputs: ExplorerInputs
    highlight_input: str | None
    highlight_change: str
    expected_result: str
    expected_route: str
    explanation: str


@dataclass(frozen=True)
class Lesson:
    """One guided question with one or more related input states."""

    key: str
    title: str
    question: str
    states: tuple[LessonState, ...]


_FULL_LAUNCH = -sqrt(7.56)


def _inputs(
    m: float, q: float, w: float, x_0: float = _FULL_LAUNCH
) -> ExplorerInputs:
    """Build exact edge speeds for the common shell and oblique ray."""
    return ExplorerInputs(1.0, 3.0, m + q, 3.0 * m + q, 1.2, x_0, w)


LESSONS: tuple[Lesson, ...] = (
    Lesson(
        "homology",
        "1. Homology reference",
        "Why is there only one potential encounter?",
        (
            LessonState(
                "initial",
                "Homologous expansion",
                _inputs(1.0, 0.0, 0.3),
                "q",
                "Reference: q = 0; radial edge velocities are 1 and 3.",
                "v∥ = x and G = 1: one forward crossing at x = 0.3, with no frequency reversal. Two equal algebraic slots represent that one encounter.",
                "zero-intercept Q=0",
                "The projected law is linear because the radial intercept vanishes. The exact Q=0 solver route repeats its root algebraically.",
            ),
        ),
    ),
    Lesson(
        "linked",
        "2. Linked nonhomology",
        "Where did three potential encounters come from?",
        (
            LessonState(
                "initial",
                "Three forward crossings",
                _inputs(1.0, -2.34375, 0.3),
                "q",
                "Change only q from 0 to −2.34375; radial edge velocities become −1.34375 and 0.65625.",
                "Two G = 0 reversals split the full ray into three monotone pieces. Three original-equation crossings are forward; the fourth finite quartic slot has the wrong unsquared sign.",
                "general quartic",
                "Squaring (m x − w)²(x² + b²) − q²x² = 0 creates four algebraic slots, but loses the sign condition from m x + q x/√(x² + b²) = w.",
            ),
        ),
    ),
    Lesson(
        "launch",
        "3. Move the launch",
        "Why do only two, then one, crossings remain?",
        (
            LessonState(
                "two_crossings",
                "Launch at x₀ = −0.5",
                _inputs(1.0, -2.34375, 0.3, -0.5),
                "x_0",
                "Change only launch x₀ to −0.5; the law and target stay fixed.",
                "The first full-ray root is behind the start. One active G = 0 turn and two forward crossings remain.",
                "general quartic",
                "The quartic slots stay tied to the same full-ray equation; the packet's forward starting point removes earlier solutions.",
            ),
            LessonState(
                "one_crossing",
                "Launch at x₀ = 1.2",
                _inputs(1.0, -2.34375, 0.3, 1.2),
                "x_0",
                "Change only launch x₀ to 1.2; the law and target stay fixed.",
                "The active outbound leg is monotone and contains one forward crossing.",
                "general quartic",
                "The two earlier physical roots are now behind the packet, even though they remain in the four solver slots.",
            ),
        ),
    ),
    Lesson(
        "target",
        "4. Move the line target",
        "When do roots disappear or meet?",
        (
            LessonState(
                "above_range",
                "Target above active range",
                _inputs(1.0, -2.34375, 0.7),
                "w",
                "Change only w from 0.3 to 0.7.",
                "The target is above the active projected-speed maximum, so there is no forward unsquared crossing.",
                "general quartic",
                "An algebraic root can still occur beyond the shell; the active projected curve never reaches this target.",
            ),
            LessonState(
                "touch",
                "Target touches the crest",
                _inputs(1.0, -2.34375, 0.50625),
                "w",
                "Change only w from 0.3 to 0.50625.",
                "The target touches the inbound crest analytically at x = −0.9 and crosses again outbound. The production solver's two nearby slots at the contact are reported separately from that analytic event.",
                "general quartic",
                "At a tangent root, a sign-change bracket alone misses the contact; compare the analytic G = 0 position with the production slots.",
            ),
        ),
    ),
    Lesson(
        "collapse",
        "5. Whole-shell collapse",
        "Can inflow still have repeated forward crossings?",
        (
            LessonState(
                "initial",
                "Negative radial speed throughout",
                _inputs(1.0, -3.5, 0.6),
                None,
                "Set m = 1, q = −3.5, w = 0.6; radial edge velocities are −2.5 and −0.5.",
                "The shell is wholly collapsing, yet the inbound ray has two independently bracketed crossings near x = −2.57170109 and −0.33120996. The remaining finite slots are outside the shell near 3.94878442 and wrong-sign near 0.15412663.",
                "general quartic",
                "Inflow describes radial motion; the projected speed along an oblique packet ray can turn twice even while radial velocity stays negative.",
            ),
        ),
    ),
    Lesson(
        "homologous_collapse",
        "6. Homologous collapse reference",
        "What does the limiting inflow do?",
        (
            LessonState(
                "initial",
                "Linear projected inflow",
                _inputs(-1.0, 0.0, 0.6),
                None,
                "Set m = −1, q = 0, w = 0.6; radial edge velocities are −1 and −3.",
                "v∥ = −x and G = −1: one forward crossing at x = −0.6 and the opposite comoving-frequency direction from homologous expansion.",
                "zero-intercept Q=0",
                "This is a collapse reference, not a one-input comparison with the preceding nonhomologous collapse case.",
            ),
        ),
    ),
    Lesson(
        "branches",
        "7. Branch and event sampler",
        "Which special routes and events change the decision path?",
        (
            LessonState(
                "radial",
                "Radial inward: H² = 0",
                PRESETS["radial inward"].inputs,
                None,
                "Choose a radial ray with b = 0.",
                "The zero impact parameter selects the radial structured-solver route.",
                "radial H²=0",
                "A radial ray has H² = 0, so the structured solver takes its radial algebraic branch before Ferrari.",
            ),
            LessonState(
                "uniform",
                "Uniform radial speed: Δv = 0",
                PRESETS["uniform speed analytic"].inputs,
                None,
                "Set equal radial speeds at the two shell edges.",
                "The constant-speed projection has one analytic oblique root.",
                "Δv=0 analytic",
                "The outer adapter handles a constant radial speed analytically before calling the structured quartic.",
            ),
            LessonState(
                "line_center",
                "Line center: N = 0",
                PRESETS["line center"].inputs,
                None,
                "Set the line target w = 0.",
                "The zero target selects the line-center structured-solver route.",
                "line-center N=0",
                "N = w/Δv is exactly zero, so the structured solver uses its line-center branch before general Ferrari factorization.",
            ),
            LessonState(
                "scatter",
                "Earlier scatter clips the flight",
                PRESETS["scatter truncation"].inputs,
                None,
                "Assume an earlier electron scatter at x = 1.6.",
                "A later algebraic crossing is excluded by event ordering; only the earlier eligible potential resonance remains.",
                "general quartic",
                "The solver still returns the full-ray roots; packet event ownership restricts the active interval.",
            ),
            LessonState(
                "positive_discriminant",
                "Positive reconstructed discriminant",
                _inputs(1.0, 0.5, 0.3),
                None,
                "Set m = 1, q = 0.5, w = 0.3 on the common ray.",
                "The active projected curve is monotone with one physical crossing, while the reconstructed Ferrari resolvent discriminant is positive.",
                "general quartic",
                "The discriminant selects a numerical factorization branch. Its sign does not directly count physical crossings; unsquared and forward-flight filters still decide that.",
            ),
        ),
    ),
)

LESSONS_BY_KEY: dict[str, Lesson] = {lesson.key: lesson for lesson in LESSONS}


def get_lesson(key: str) -> Lesson:
    """Return the guided lesson with the specified stable key.

    Parameters
    ----------
    key : str
        Stable lesson key.

    Returns
    -------
    Lesson
        Matching lesson.
    """
    return LESSONS_BY_KEY[key]


def get_lesson_state(lesson_key: str, state_key: str) -> LessonState:
    """Return one exact-input state within a guided lesson.

    Parameters
    ----------
    lesson_key : str
        Stable lesson key.
    state_key : str
        Stable state key inside the lesson.

    Returns
    -------
    LessonState
        Matching state.
    """
    return next(
        state
        for state in get_lesson(lesson_key).states
        if state.key == state_key
    )
