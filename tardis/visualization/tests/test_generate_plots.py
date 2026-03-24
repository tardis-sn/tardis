"""Tests for the generate_plots plotting script.

These are *unit* tests: the TARDIS simulation is replaced with a lightweight
mock object so the tests run quickly and without requiring any data files.
Integration tests against a real simulation are left to CI with the usual
``simulation_simple_tracked`` fixture.
"""

import os
import sys
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Add the docs directory to sys.path so generate_plots can be imported
# without installing it as a package.
# ---------------------------------------------------------------------------
_DOCS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../../../docs/analyzing_tardis/visualization")
)
if _DOCS_DIR not in sys.path:
    sys.path.insert(0, _DOCS_DIR)

import generate_plots as _mod  # noqa: E402 – intentional late import after path fix


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_mock_sim():
    """Return a mock Simulation with the minimal attributes needed by generate_plots."""
    spectrum = MagicMock()
    spectrum.wavelength = [5000.0, 6000.0, 7000.0]
    spectrum.luminosity_density_lambda = [1.0e10, 2.0e10, 1.5e10]

    sim = MagicMock()
    sim.spectrum_solver.spectrum = spectrum
    return sim


@pytest.fixture()
def mock_sim():
    return _make_mock_sim()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestGeneratePlotsScript:
    """Unit tests for generate_plots.py."""

    def test_module_exposes_generate_plots(self):
        """The module must expose a public callable named generate_plots."""
        assert callable(_mod.generate_plots)

    def test_output_directory_created_when_missing(self, tmp_path, mock_sim):
        """generate_plots must create the output directory if it does not exist."""
        output_dir = str(tmp_path / "new_dir")
        assert not os.path.exists(output_dir)

        with patch("generate_plots.run_tardis", return_value=mock_sim), \
             patch("generate_plots.SDECPlotter") as mock_sdec, \
             patch("generate_plots.tardis") as mock_tardis:

            mock_tardis.__version__ = "0.0.0"

            ax = MagicMock()
            mock_sdec.from_simulation.return_value.generate_plot_mpl.return_value = ax

            # Prevent LIVPlotter import from reaching real tardis
            with patch.dict("sys.modules", {"tardis.visualization": MagicMock(spec=[])}):
                _mod.generate_plots("dummy.yml", output_dir=output_dir)

        assert os.path.isdir(output_dir)

    def test_sdec_plot_saved_with_correct_filename(self, tmp_path, mock_sim):
        """generate_plots should call savefig with the expected file path."""
        output_dir = str(tmp_path / "plots")

        with patch("generate_plots.run_tardis", return_value=mock_sim), \
             patch("generate_plots.SDECPlotter") as mock_sdec, \
             patch("generate_plots.tardis") as mock_tardis:

            mock_tardis.__version__ = "0.0.0"
            ax = MagicMock()
            mock_sdec.from_simulation.return_value.generate_plot_mpl.return_value = ax

            with patch.dict("sys.modules", {"tardis.visualization": MagicMock(spec=[])}):
                _mod.generate_plots(
                    "dummy.yml",
                    output_dir=output_dir,
                    prefix="run01",
                    output_format="png",
                )

        expected_path = os.path.join(output_dir, "run01_sdec.png")
        ax.figure.savefig.assert_called_once_with(expected_path, bbox_inches="tight")

    def test_fallback_spectrum_saved_when_sdec_fails(self, tmp_path, mock_sim):
        """When SDECPlotter raises, the fallback spectrum file must be saved."""
        output_dir = str(tmp_path / "fallback")

        with patch("generate_plots.run_tardis", return_value=mock_sim), \
             patch("generate_plots.SDECPlotter", side_effect=RuntimeError("mock failure")), \
             patch("generate_plots.tardis") as mock_tardis, \
             patch("generate_plots.plt") as mock_plt:

            mock_tardis.__version__ = "0.0.0"
            fig_mock = MagicMock()
            mock_plt.subplots.return_value = (fig_mock, MagicMock())

            with patch.dict("sys.modules", {"tardis.visualization": MagicMock(spec=[])}):
                _mod.generate_plots("dummy.yml", output_dir=output_dir, prefix="fb")

        fallback_path = os.path.join(output_dir, "fb_spectrum_fallback.pdf")
        fig_mock.savefig.assert_called_once_with(fallback_path, bbox_inches="tight")
