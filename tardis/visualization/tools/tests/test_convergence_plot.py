"""Tests for Convergence Plots."""

from copy import deepcopy

import pytest
from tardis.tests.test_util import monkeysession
from tardis import run_tardis
from tardis.visualization.tools.convergence_plot import (
    ConvergencePlots,
    transition_colors,
)
from collections import defaultdict
import plotly.graph_objects as go
from astropy import units as u


@pytest.fixture(scope="module", params=[0, 1, 2])
def convergence_plots(request):
    """Initialize ConvergencePlots class and build empty plots."""
    convergence_plots = ConvergencePlots(iterations=request.param)
    convergence_plots.build(display_plot=False)
    return convergence_plots


@pytest.fixture()
def fetch_luminosity_data(convergence_plots):
    """Prepare data for t_inner and luminosity plot."""
    for item in [2] * convergence_plots.iterations:
        convergence_plots.fetch_data(
            name="t_inner", value=item, item_type="value"
        )
        for item2 in convergence_plots.luminosities:
            convergence_plots.fetch_data(
                name=item2, value=item, item_type="value"
            )


def test_transition_colors():
    """Test whether the object returned by the transition_colors function is a list of appropriate length."""
    iterations = 3
    colors = transition_colors(length=iterations)
    assert type(colors) == list
    assert len(colors) == iterations


def test_convergence_construction(convergence_plots):
    """Test the construction of the ConvergencePlots class."""
    assert convergence_plots.iterable_data == {}
    assert convergence_plots.value_data == defaultdict(list)
    assert convergence_plots.luminosities == [
        "Emitted",
        "Absorbed",
        "Requested",
    ]


def test_fetch_data(convergence_plots):
    """Test values of variables updated by fetch_data function."""
    convergence_plots.fetch_data(
        name="iterable", value=range(3), item_type="iterable"
    )
    convergence_plots.fetch_data(name="value", value=0, item_type="value")

    assert convergence_plots.iterable_data["iterable"] == range(3)
    assert convergence_plots.value_data["value"] == [0]


def test_build(convergence_plots):
    """Test if convergence plots are instances of plotly.graph_objs.FigureWidget() and have appropriate number of traces."""
    assert type(convergence_plots.plasma_plot) == go.FigureWidget
    assert type(convergence_plots.t_inner_luminosities_plot) == go.FigureWidget

    # check number of traces
    assert len(convergence_plots.t_inner_luminosities_plot.data) == 5
    assert len(convergence_plots.plasma_plot.data) == 2


@pytest.mark.usefixtures("fetch_luminosity_data")
def test_update_t_inner_luminosities_plot(convergence_plots):
    """Test the number of traces and length of x and y values."""
    n_iterations = convergence_plots.iterations
    convergence_plots.update_t_inner_luminosities_plot()

    # check number of traces
    assert len(convergence_plots.t_inner_luminosities_plot.data) == 5

    for index in range(0, 5):
        # check x and y values for all traces
        assert (
            len(convergence_plots.t_inner_luminosities_plot.data[index].x)
            == n_iterations
        )
        assert (
            len(convergence_plots.t_inner_luminosities_plot.data[index].y)
            == n_iterations
        )

    # check range of x axes
    for axis in ["xaxis", "xaxis2", "xaxis3"]:
        convergence_plots.t_inner_luminosities_plot["layout"][axis]["range"] = [
            0,
            convergence_plots.iterations + 1,
        ]


@pytest.mark.usefixtures("fetch_luminosity_data")
def test_update_plasma_plots(convergence_plots):
    """Test the state of plasma plots after updating."""
    n_iterations = convergence_plots.iterations
    expected_n_traces = 2 * n_iterations + 2
    velocity = range(0, n_iterations) * u.m / u.s

    convergence_plots.fetch_data(
        name="velocity", value=velocity, item_type="iterable"
    )

    w_val = list(range(n_iterations))
    t_rad_val = [item * 2 for item in w_val]

    for _ in range(n_iterations):
        convergence_plots.fetch_data(
            name="t_rad",
            value=t_rad_val,
            item_type="iterable",
        )
        convergence_plots.fetch_data(
            name="w",
            value=w_val,
            item_type="iterable",
        )
        convergence_plots.update_plasma_plots()
        convergence_plots.current_iteration += 1

    # check number of traces
    assert len(convergence_plots.plasma_plot.data) == expected_n_traces

    # traces are added alternatively
    # trace 0 and 1 and empty
    assert convergence_plots.plasma_plot.data[0].x == None
    assert convergence_plots.plasma_plot.data[1].x == None

    # check other traces
    for index in list(range(expected_n_traces))[::2][1:]:
        # check values for t_rad subplot
        assert convergence_plots.plasma_plot.data[index].xaxis == "x"
        assert convergence_plots.plasma_plot.data[index].yaxis == "y"
        assert (
            convergence_plots.plasma_plot.data[index].y[:-1] == tuple(t_rad_val)
        ).all()
        assert convergence_plots.plasma_plot.data[index].x == tuple(
            velocity.to(u.km / u.s).value
        )

    for index in list(range(expected_n_traces))[1::2][1:]:
        # check values for w subplot
        assert convergence_plots.plasma_plot.data[index].xaxis == "x2"
        assert convergence_plots.plasma_plot.data[index].yaxis == "y2"
        assert (
            convergence_plots.plasma_plot.data[index].y[:-1] == tuple(w_val)
        ).all()
        assert convergence_plots.plasma_plot.data[index].x == tuple(
            velocity.to(u.km / u.s).value
        )


@pytest.mark.usefixtures("fetch_luminosity_data")
def test_override_plot_parameters(convergence_plots):
    """Test if default plot properties are overridden properly."""
    parameters = {
        "data": {
            "line": {"dash": "dot"},
            "mode": "lines+markers",
        },
        "layout": {"xaxis2": {"showgrid": False}},
    }
    convergence_plots.override_plot_parameters(
        convergence_plots.plasma_plot, parameters=parameters
    )
    convergence_plots.override_plot_parameters(
        convergence_plots.t_inner_luminosities_plot, parameters=parameters
    )

    # data properties will be applied across all traces equally
    # testing plot parameters of t_inner and luminosity plot
    for trace_index in range(5):
        assert (
            convergence_plots.t_inner_luminosities_plot.data[
                trace_index
            ].line.dash
            == "dot"
        )
        assert (
            convergence_plots.t_inner_luminosities_plot.data[trace_index].mode
            == "lines+markers"
        )
    # checking layout
    assert (
        convergence_plots.t_inner_luminosities_plot["layout"]["xaxis2"][
            "showgrid"
        ]
        == False
    )

    # testing plot parameters for plasma plot
    for trace_index in range(2):
        assert (
            convergence_plots.plasma_plot.data[trace_index].line.dash == "dot"
        )
        assert (
            convergence_plots.plasma_plot.data[trace_index].mode
            == "lines+markers"
        )
    # checking layout for plasma plot
    assert (
        convergence_plots.plasma_plot["layout"]["xaxis2"]["showgrid"] == False
    )


def test_convergence_plot_command_line(
    config_verysimple, atomic_dataset, monkeysession
):
    monkeysession.setattr(
        "tardis.simulation.base.is_notebook",
        lambda: False,
    )
    atomic_data = deepcopy(atomic_dataset)
    with pytest.raises(RuntimeError):
        run_tardis(
            config_verysimple,
            atom_data=atomic_data,
            show_convergence_plots=True,
        )
