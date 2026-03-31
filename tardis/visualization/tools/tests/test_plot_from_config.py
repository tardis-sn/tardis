"""Regression tests for plot generation from config."""

import os
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from tardis.visualization.tools.plot_from_config import (
    build_arg_parser,
    save_plots,
)
from tardis.visualization.tools.sdec_plot import SDECPlotter


@pytest.fixture
def simulation_for_plots(simulation_simple_tracked):
    """Provide a tracked simulation for plotting tests."""
    return simulation_simple_tracked


def test_build_arg_parser_defaults():
    """Test parser defaults."""
    parser = build_arg_parser()
    args = parser.parse_args(["config.yml"])

    assert args.config == "config.yml"
    assert args.atom_data is None
    assert args.output_prefix is None
    assert args.output_dir == "."
    assert args.format == "png"


def test_save_plots_png_regression(simulation_for_plots, regression_data):
    """Generate PNG files and persist deterministic metadata as regression data."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        output_prefix = "plot_from_config"
        save_plots(
            simulation_for_plots,
            [SDECPlotter],
            output_prefix=output_prefix,
            output_dir=tmp_dir,
            fmt="png",
        )

        saved_files = sorted(str(path) for path in Path(tmp_dir).glob("*.png"))

        metadata = pd.DataFrame(
            {
                "basename": sorted(
                    os.path.basename(path) for path in saved_files
                ),
                "exists": [os.path.exists(path) for path in sorted(saved_files)],
                "non_empty": [
                    os.path.getsize(path) > 0 for path in sorted(saved_files)
                ],
            }
        )

        expected = regression_data.sync_dataframe(
            metadata, key="save_plots_png_metadata"
        )
        pd.testing.assert_frame_equal(metadata, expected)


def test_save_plots_pdf_regression(simulation_for_plots, regression_data):
    """Generate PDF files and persist deterministic metadata as regression data."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        output_prefix = "plot_from_config"
        save_plots(
            simulation_for_plots,
            [SDECPlotter],
            output_prefix=output_prefix,
            output_dir=tmp_dir,
            fmt="pdf",
        )

        saved_files = sorted(str(path) for path in Path(tmp_dir).glob("*.pdf"))

        metadata = pd.DataFrame(
            {
                "basename": sorted(
                    os.path.basename(path) for path in saved_files
                ),
                "exists": [os.path.exists(path) for path in sorted(saved_files)],
                "non_empty": [
                    os.path.getsize(path) > 0 for path in sorted(saved_files)
                ],
            }
        )

        expected = regression_data.sync_dataframe(
            metadata, key="save_plots_pdf_metadata"
        )
        pd.testing.assert_frame_equal(metadata, expected)
