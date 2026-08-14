"""Tests for publication-ready atomic-data summary tables."""

import json

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from tardis.io.atom_data import AtomData, export_atom_data_summary


# Test 1: Can every built in journal format generate the tex output?
@pytest.mark.parametrize("journal", ["aa", "aas", "mnras", "nature"])
def test_export_atom_data_summary_supports_builtin_journals(
    tmp_path,
    journal,
):
    # Create the smallest atomic dataset that can produce a table row.
    levels = pd.DataFrame({"atomic_number": [1, 1], "ion_number": [0, 0]})
    lines = pd.DataFrame({"atomic_number": [1], "ion_number": [0]})
    atom_data = object.__new__(AtomData)
    atom_data.levels = levels
    atom_data.lines = lines

    # Export the same data with each built-in journal style.
    export_atom_data_summary(
        atom_data,
        tmp_path / journal,
        journal=journal,
    )

    # Check that each style creates a LaTeX file.
    assert (tmp_path / f"{journal}.tex").is_file()

# Test 2: Are ion ranges, level counts, and totals in the table accurate?
def test_export_atom_data_summary_from_atom_data(tmp_path):
    # Include two usable hydrogen ions and one unusable helium ion.
    levels = pd.DataFrame(
        {
            "atomic_number": [1, 1, 1, 1, 1, 2],
            "ion_number": [0, 0, 1, 1, 3, 0],
        }
    )
    lines = pd.DataFrame(
        {
            "atomic_number": [1, 1, 1, 2],
            "ion_number": [0, 0, 1, 0],
        }
    )
    atom_data = object.__new__(AtomData)
    atom_data.levels = levels
    atom_data.lines = lines

    # Generate the summary directly from an AtomData object.
    output_path = tmp_path / "summary"
    summary = export_atom_data_summary(
        atom_data,
        output_path,
        journal="mnras",
    )

    # Check the counts, ion range, total, and written files.
    expected = pd.DataFrame(
        {
            "Element": [r"$\mathrm{H}_{1}$", "Total"],
            "Ion stages": ["I--II", ""],
            "No. of Levels": [4, 4],
            "No. of Lines": [3, 3],
        }
    )
    assert_frame_equal(summary, expected)
    assert output_path.with_suffix(".txt").is_file()
    latex = output_path.with_suffix(".tex").read_text(encoding="utf-8")
    assert r"\documentclass{article}" in latex
    assert r"$\mathrm{H}_{1}$ & I--II & 4 & 3 \\" in latex
    assert r"\textbf{Total}" in latex

# Test 3: Do additional features (element filtering and custom journal config) work as expected?
def test_export_atom_data_summary_with_elements_only(tmp_path):
    # Create usable data for hydrogen and helium.
    levels = pd.DataFrame(
        {
            "atomic_number": [1, 1, 2, 2],
            "ion_number": [0, 0, 0, 0],
        }
    )
    lines = pd.DataFrame(
        {
            "atomic_number": [1, 2],
            "ion_number": [0, 0],
        }
    )
    atom_data = object.__new__(AtomData)
    atom_data.levels = levels
    atom_data.lines = lines

    # Export only hydrogen using the default journal.
    output_path = tmp_path / "hydrogen_table"
    summary = export_atom_data_summary(
        atom_data,
        output_path,
        elements=["H"],
    )

    # Check that only hydrogen rows are present.
    assert summary["Element"].tolist() == [r"$\mathrm{H}_{1}$", "Total"]


def test_export_atom_data_summary_with_custom_journal(tmp_path):
    # Create usable data for hydrogen and helium.
    levels = pd.DataFrame(
        {
            "atomic_number": [1, 1, 2, 2],
            "ion_number": [0, 0, 0, 0],
        }
    )
    lines = pd.DataFrame(
        {
            "atomic_number": [1, 2],
            "ion_number": [0, 0],
        }
    )
    atom_data = object.__new__(AtomData)
    atom_data.levels = levels
    atom_data.lines = lines
    # Define a minimal custom journal template without a bold total.
    config_path = tmp_path / "journal.json"
    config_path.write_text(
        json.dumps(
            {
                "bold_total": False,
                "template": "BEGIN\n{{TABLE_ROWS}}\nEND\n",
            }
        ),
        encoding="utf-8",
    )

    # Export using the custom template.
    output_path = tmp_path / "custom_table.tex"
    summary = export_atom_data_summary(
        atom_data,
        output_path,
        journal=config_path,
    )

    # Check custom formatting is applied and total is not bold.
    latex = output_path.with_suffix(".tex").read_text(encoding="utf-8")
    assert latex.startswith("BEGIN\n")
    assert r"\textbf{Total}" not in latex

# Test 4: Can the original HDF files be read directly, not just the AtomData object?-spli
def test_export_atom_data_summary_from_hdf(tmp_path):
    # Store a minimal silicon dataset in the expected HDF tables.
    input_path = tmp_path / "atom_data.h5"
    pd.DataFrame(
        {
            "atomic_number": [14, 14],
            "ion_number": [1, 1],
        }
    ).to_hdf(input_path, key="levels_data")
    pd.DataFrame(
        {
            "atomic_number": [14],
            "ion_number": [1],
        }
    ).to_hdf(input_path, key="lines_data")

    # Generate a summary directly from the HDF file.
    summary = export_atom_data_summary(
        input_path,
        tmp_path / "silicon",
        journal="aas",
    )

    # Check that the HDF data is counted and labeled correctly.
    assert summary.iloc[0].to_dict() == {
        "Element": r"$\mathrm{Si}_{14}$",
        "Ion stages": "II",
        "No. of Levels": 2,
        "No. of Lines": 1,
    }
