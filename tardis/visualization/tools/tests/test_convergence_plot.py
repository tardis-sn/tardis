from attr import s
import pytest
from tardis.visualization.tools.convergence_plot import (
    ConvergencePlots,
    transition_colors,
)
from collections import defaultdict
import plotly.graph_objects as go
from astropy import units as u


@pytest.fixture(scope="module", params=[0, 1, 2])
def convergence_class(request):
    """
    Fixture to initialize ConvergencePlots class and build empty plots.

    Returns
    ----------
    convergence_class : ConvergencePlots
    """
    convergence_plots = ConvergencePlots(iterations=request.param)
    convergence_plots.build(display_plot=False)
    return convergence_plots


@pytest.fixture()
def fetch_luminosity_data(convergence_class):
    """
    Prepares data for luminosity plot.

    Parameters
    ----------
    ConvergencePlots class

    """
    for item in [2] * convergence_class.iterations:
        convergence_class.fetch_data(
            name="t_inner", value=item, item_type="value"
        )
        for item2 in convergence_class.luminosities:
            convergence_class.fetch_data(
                name=item2, value=item, item_type="value"
            )


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
    Test if convergence plots are instances of plotly.graph_objs.FigureWidget() and
    have appropriate number of traces.

    Parameters
    ----------
    convergence_class : ConvergencePlots
    """
    assert type(convergence_class.plasma_plot) == go.FigureWidget
    assert type(convergence_class.luminosity_plot) == go.FigureWidget

    # check number of traces
    assert len(convergence_class.luminosity_plot.data) == 5
    assert len(convergence_class.plasma_plot.data) == 2


@pytest.mark.usefixtures("fetch_luminosity_data")
def test_update_luminosity_plot(convergence_class):
    """
    Tests the number of traces and length of x and y values.

    Parameters
    ----------
    convergence_class : ConvergencePlots
    """
    n_iterations = convergence_class.iterations
    convergence_class.update_luminosity_plot()

    # check number of traces
    assert len(convergence_class.luminosity_plot.data) == 5

    for index in range(0, 5):
        # check x and y values for all traces
        assert (
            len(convergence_class.luminosity_plot.data[index].x) == n_iterations
        )
        assert (
            len(convergence_class.luminosity_plot.data[index].y) == n_iterations
        )
    for axis in ["xaxis", "xaxis2", "xaxis3"]:
        convergence_class.luminosity_plot["layout"][axis]["range"] = [
            0,
            convergence_class.iterations + 1,
        ]


@pytest.mark.usefixtures("fetch_luminosity_data")
def test_update_plasma_plots(convergence_class):
    """
    Tests the state of plasma plots after updating.

    Parameters
    ----------
    convergence_class : ConvergencePlots
    """

    n_iterations = convergence_class.iterations
    expected_n_traces = 2 * n_iterations + 2
    velocity = range(0, n_iterations) * u.m / u.s

    convergence_class.fetch_data(
        name="velocity", value=velocity, item_type="iterable"
    )
    w_val = list(range(n_iterations))
    t_rad_val = [item * 2 for item in w_val]

    for _ in range(n_iterations):
        convergence_class.fetch_data(
            name="t_rad",
            value=t_rad_val,
            item_type="iterable",
        )
        convergence_class.fetch_data(
            name="w",
            value=w_val,
            item_type="iterable",
        )
        convergence_class.update_plasma_plots()
        convergence_class.current_iteration += 1

    # check number of traces
    assert len(convergence_class.plasma_plot.data) == expected_n_traces

    # check if the first two traces are empty
    assert convergence_class.plasma_plot.data[0].x == None
    assert convergence_class.plasma_plot.data[1].x == None

    for index in list(range(expected_n_traces))[::2][1:]:
        # check values for t_rad subplot
        assert convergence_class.plasma_plot.data[index].xaxis == "x"
        assert convergence_class.plasma_plot.data[index].yaxis == "y"
        assert convergence_class.plasma_plot.data[index].y == tuple(t_rad_val)
        assert convergence_class.plasma_plot.data[index].x == tuple(
            velocity.to(u.km / u.s).value
        )

    for index in list(range(expected_n_traces))[1::2][1:]:
        # check values for w subplot
        assert convergence_class.plasma_plot.data[index].xaxis == "x2"
        assert convergence_class.plasma_plot.data[index].yaxis == "y2"
        assert convergence_class.plasma_plot.data[index].y == tuple(w_val)
        assert convergence_class.plasma_plot.data[index].x == tuple(
            velocity.to(u.km / u.s).value
        )


@pytest.mark.usefixtures("fetch_luminosity_data")
def test_override_plot_parameters(convergence_class):
    parameters = {
        "data": {
            "line": {"dash": "dot"},
            "mode": "lines+markers",
        },
        "layout": {"xaxis2": {"showgrid": False}},
    }
    convergence_class.override_plot_parameters(
        convergence_class.plasma_plot, parameters=parameters
    )
    convergence_class.override_plot_parameters(
        convergence_class.luminosity_plot, parameters=parameters
    )

    # testing plot parameters of luminosity plot
    for trace_index in range(5):
        assert (
            convergence_class.luminosity_plot.data[trace_index].line.dash
            == "dot"
        )
        assert (
            convergence_class.luminosity_plot.data[trace_index].mode
            == "lines+markers"
        )
    assert (
        convergence_class.luminosity_plot["layout"]["xaxis2"]["showgrid"]
        == False
    )

    # testing plot parameters for plasma plot
    for trace_index in range(2):
        assert (
            convergence_class.plasma_plot.data[trace_index].line.dash == "dot"
        )
        assert (
            convergence_class.plasma_plot.data[trace_index].mode
            == "lines+markers"
        )
    assert (
        convergence_class.plasma_plot["layout"]["xaxis2"]["showgrid"] == False
    )
