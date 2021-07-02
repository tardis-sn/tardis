import pytest
from tardis.visualization.tools.convergence_plot import (
    ConvergencePlots,
    transition_colors,
)
from collections import defaultdict
import plotly.graph_objects as go


@pytest.fixture
def convergence_class() -> ConvergencePlots:
    return ConvergencePlots(iterations=3)


def test_transition_colors():
    """
    Test whether the object returned by the
    transition_colors function is a list of appropriate length.
    """
    iterations = 3
    colors = transition_colors(length=iterations)
    assert type(colors) == list
    assert len(colors) == iterations


def test_convergence_construction(convergence_class):
    """
    Test the construction of the ConvergencePlots class.

    Parameters
    ----------
    convergence_class : ConvergencePlots
    """

    assert convergence_class.iterable_data == {}
    assert convergence_class.value_data == defaultdict(list)
    assert convergence_class.luminosities == [
        "Emitted",
        "Absorbed",
        "Requested",
    ]


def test_fetch_data(convergence_class):
    """
    Test values of variables updated by fetch_data function.

    Parameters
    ----------
    convergence_class : ConvergencePlots
    """
    convergence_class.fetch_data(
        name="iterable", value=range(3), item_type="iterable"
    )
    convergence_class.fetch_data(name="value", value=0, item_type="value")

    assert convergence_class.iterable_data["iterable"] == range(3)
    assert convergence_class.value_data["value"] == [0]


def test_build(convergence_class):
    """
    Test if convergence plots are instances of plotly.graph_objs.FigureWidget()

    Parameters
    ----------
    convergence_class : ConvergencePlots
    """
    convergence_class.build(display_plot=False)
    assert type(convergence_class.plasma_plot) == go.FigureWidget
    assert type(convergence_class.luminosity_plot) == go.FigureWidget

    # check number of traces
    assert len(convergence_class.luminosity_plot.data) == 5
    assert len(convergence_class.plasma_plot.data) == 2
