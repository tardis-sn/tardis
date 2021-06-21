import pytest
from tardis.visualization.tools.convergence_plot import ConvergencePlots, transistion_colors
from collections import defaultdict


@pytest.fixture
def convergence_class() -> ConvergencePlots:
    return ConvergencePlots(iterations=3)


def test_transition_colors():
    iterations = 3
    colors = transistion_colors(iterations=iterations)
    assert type(colors) == list
    assert len(colors) == iterations


def test_convergence_construction(convergence_class):
    assert convergence_class.iterable_data == {}
    assert convergence_class.value_data == defaultdict(list)
    assert convergence_class.luminosities == ["Emitted", "Absorbed", "Requested"]


def test_fetch_data(convergence_class):
    convergence_class.fetch_data(name="iterable", value=range(3), type="iterable")
    convergence_class.fetch_data(name="value", value=0, type="value")
    assert len(convergence_class.iterable_data["iterable"]) == 3
