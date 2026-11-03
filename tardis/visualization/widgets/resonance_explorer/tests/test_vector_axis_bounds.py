"""Fixed physical coordinates in the guided vector construction."""

import numpy.testing as npt
from matplotlib.patches import FancyArrowPatch

from tardis.visualization.widgets.resonance_explorer.analysis import (
    analyze_case,
)
from tardis.visualization.widgets.resonance_explorer.lessons import (
    get_lesson_state,
)
from tardis.visualization.widgets.resonance_explorer.render import (
    render_guided_case,
)


def test_guided_vector_coordinates_match_shell_for_each_launch() -> None:
    """Claim: shortening a packet flight does not rescale its vector diagram.

    Regime: the linked velocity law and target with two and one forward
    crossings as launch x0 moves. Verification: both views use the analytic
    shell extent +/-1.08*r_out in x and y, independent of the active ray.
    This detects bounds inferred from visible vector or trajectory endpoints.
    """
    for state_key in ("two_crossings", "one_crossing"):
        inputs = get_lesson_state("launch", state_key).inputs
        figure = render_guided_case(analyze_case(inputs))
        figure.canvas.draw()
        shell = next(
            axis
            for axis in figure.axes
            if axis.get_gid() == "guided shell geometry"
        )
        vectors = next(
            axis
            for axis in figure.axes
            if axis.get_gid() == "guided vector construction"
        )
        half_width = 1.08 * inputs.r_out
        expected_limits = (-half_width, half_width)
        npt.assert_allclose(shell.get_xlim(), expected_limits, rtol=0, atol=0)
        npt.assert_allclose(shell.get_ylim(), expected_limits, rtol=0, atol=0)
        npt.assert_allclose(vectors.get_xlim(), expected_limits, rtol=0, atol=0)
        npt.assert_allclose(vectors.get_ylim(), expected_limits, rtol=0, atol=0)
        npt.assert_allclose(
            vectors.get_window_extent().width,
            shell.get_window_extent().width,
            rtol=0,
            atol=1e-9,
        )
        npt.assert_allclose(
            vectors.get_window_extent().height,
            shell.get_window_extent().height,
            rtol=0,
            atol=1e-9,
        )
        npt.assert_allclose(
            vectors.transData.transform((0, 0))[0],
            shell.transData.transform((0, 0))[0],
            rtol=0,
            atol=1e-9,
        )


def test_fixed_guided_frame_contains_velocity_arrow_tips() -> None:
    """Claim: the fixed shell frame still contains the constructed vectors.

    Regime: an outbound crossing and a target tangent near the outer edge.
    Verification: each mathematical arrow tip lies inside the analytic shell
    frame. This catches display scales fitted to a changing flight extent.
    """
    for lesson_key, state_key in (
        ("launch", "one_crossing"),
        ("target", "touch"),
    ):
        case = analyze_case(get_lesson_state(lesson_key, state_key).inputs)
        figure = render_guided_case(case)
        vectors = next(
            axis
            for axis in figure.axes
            if axis.get_gid() == "guided vector construction"
        )
        x_min, x_max = vectors.get_xlim()
        y_min, y_max = vectors.get_ylim()
        for artist in vectors.get_children():
            if not isinstance(artist, FancyArrowPatch):
                continue
            gid = artist.get_gid()
            if gid is None or not gid.startswith("vector-"):
                continue
            x_tip, y_tip = artist._posA_posB[1]
            assert x_min <= x_tip <= x_max, (lesson_key, state_key, gid)
            assert y_min <= y_tip <= y_max, (lesson_key, state_key, gid)
