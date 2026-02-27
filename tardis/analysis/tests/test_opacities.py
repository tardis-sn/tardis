import importlib.util
from pathlib import Path

import astropy.units as u
import pytest


MODULE_PATH = Path(__file__).resolve().parents[1] / "opacities.py"
SPEC = importlib.util.spec_from_file_location("tardis_analysis_opacities", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)
OpacityCalculator = MODULE.OpacityCalculator


class DummyModel:
    pass


@pytest.mark.parametrize("nbins", [0, -1, 1.5, True])
def test_opacity_calculator_rejects_invalid_nbins(nbins):
    with pytest.raises(ValueError, match="nbins must be a positive integer"):
        OpacityCalculator(DummyModel(), nbins=nbins)


def test_opacity_calculator_rejects_none_model():
    with pytest.raises(ValueError, match="mdl cannot be None"):
        OpacityCalculator(None)


@pytest.mark.parametrize(
    ("lam_min", "lam_max"),
    [
        (2000 * u.AA, 1000 * u.AA),
        (1000 * u.AA, 1000 * u.AA),
    ],
)
def test_opacity_calculator_rejects_invalid_wavelength_order(lam_min, lam_max):
    with pytest.raises(ValueError, match="lam_min must be smaller than lam_max"):
        OpacityCalculator(DummyModel(), lam_min=lam_min, lam_max=lam_max)


@pytest.mark.parametrize("boundary_name", ["lam_min", "lam_max"])
def test_opacity_calculator_rejects_none_wavelength_boundaries(boundary_name):
    kwargs = {boundary_name: None}

    with pytest.raises(ValueError, match=rf"{boundary_name} cannot be None"):
        OpacityCalculator(DummyModel(), **kwargs)
