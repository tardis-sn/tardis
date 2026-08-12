"""Create publication-ready summaries of TARDIS atomic data."""

from __future__ import annotations

import json
from collections.abc import Iterable
from importlib.resources import files
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from tardis.util.base import (
    atomic_number2element_symbol,
    element_symbol2atomic_number,
    int_to_roman,
)

if TYPE_CHECKING:
    from tardis.io.atom_data.base import AtomData

BUILTIN_JOURNALS = ("aa", "aas", "mnras", "nature")
JOURNAL_TEMPLATES_FILE = "atomic_data_summary_templates.json"

# Add one, de-duplicate stages, combine consecutive values into ranges, use Roman-numeral converter.
def _format_ion_stage_range(ion_numbers: Iterable[int]) -> str:
    """Format zero-based ion numbers as spectroscopic ranges."""
    # Convert the zero-based ion numbers to spectroscopic stages.
    stages = sorted({int(ion_number) + 1 for ion_number in ion_numbers})
    runs = []
    start = previous = stages[0]

    # Group consecutive stages into ranges.
    for stage in stages[1:]:
        if stage != previous + 1:
            runs.append((start, previous))
            start = stage
        previous = stage
    runs.append((start, previous))

    # Write each stage or range with Roman numerals.
    return ", ".join(
        int_to_roman(start)
        if start == end
        else f"{int_to_roman(start)}--{int_to_roman(end)}"
        for start, end in runs
    )


def _get_index_columns(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Return atomic and ion numbers as ordinary DataFrame columns."""
    # Some TARDIS tables keep these values in the index instead of columns.
    required_columns = {"atomic_number", "ion_number"}
    if required_columns.issubset(dataframe.columns):
        return dataframe
    return dataframe.reset_index()


def build_atom_data_summary(
    levels: pd.DataFrame,
    lines: pd.DataFrame,
    elements: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Summarize usable atomic levels and lines by element.

    Parameters
    ----------
    levels : pandas.DataFrame
        TARDIS atomic level data containing atomic and ion numbers.
    lines : pandas.DataFrame
        TARDIS atomic line data containing atomic and ion numbers.
    elements : collections.abc.Iterable[str], optional
        Element symbols to include. All available elements are included by
        default.

    Returns
    -------
    pandas.DataFrame
        Element, ion-stage, level-count, and line-count columns, followed by a
        total row. Only ions having more than one level and at least one line
        are included.
    """
    # TARDIS data could be stored as columns or index levels.
    # Make the atomic and ion numbers available as regular columns.
    levels = _get_index_columns(levels)
    lines = _get_index_columns(lines)

    # Combine the level and line counts for each ion.
    grouping_columns = ["atomic_number", "ion_number"]
    counts_by_ion = pd.concat(
        [
            levels.groupby(grouping_columns).size().rename("Levels"),
            lines.groupby(grouping_columns).size().rename("Lines"),
        ],
        axis=1,
    )
    # Keep ions that have usable level and line data.
    # Retain ions with more than one level and at least one line.
    eligible_ions = counts_by_ion[
        (counts_by_ion["Levels"] > 1) & counts_by_ion["Lines"].notna()
    ].reset_index()

    if elements is not None:
        # Limit the table to the requested chemical elements.
        # User-given symbols are converted to atomic numbers using existing TARDIS functionality.
        atomic_numbers = {
            element_symbol2atomic_number(element) for element in elements
        }
        eligible_ions = eligible_ions[
            eligible_ions["atomic_number"].isin(atomic_numbers)
        ]

    # Add together the ion counts for each element.
    grouped_ions = eligible_ions.groupby("atomic_number", sort=True)
    summary = pd.concat(
        [
            grouped_ions["ion_number"]
            .apply(_format_ion_stage_range)
            .rename("Ion stages"),
            grouped_ions[["Levels", "Lines"]].sum(),
        ],
        axis=1,
    )
    # Format each element as a LaTeX symbol with its atomic number.
    summary.insert(
        0,
        "Element",
        [
            rf"$\mathrm{{{atomic_number2element_symbol(atomic_number)}}}"
            rf"_{{{atomic_number}}}$"
            for atomic_number in summary.index
        ],
    )
    # Finish the table with clean integer counts and a total row.
    summary = summary.reset_index(drop=True)
    summary[["Levels", "Lines"]] = summary[["Levels", "Lines"]].astype(int)
    summary.loc[len(summary)] = [
        "Total",
        "",
        int(summary["Levels"].sum()),
        int(summary["Lines"].sum()),
    ]
    return summary


# Function accepts two input forms of HDF path or AtomData object.
def _read_atom_data_tables(
    atom_data: AtomData | str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read level and line tables from an AtomData object or HDF file."""
    if isinstance(atom_data, (str, Path)):
        # Read the original Carsus tables directly from an HDF file.
        with pd.HDFStore(atom_data, mode="r") as store:
            return store["levels_data"], store["lines_data"]

    # Use tables already loaded by TARDIS when an AtomData object is given.
    return atom_data.levels, atom_data.lines


def _load_journal_config(journal: str | Path) -> dict[str, object]:
    """Load and validate a built-in or custom journal configuration."""
    journal_path = Path(journal)

    if journal_path.suffix.lower() == ".json":
        # Load a user-provided journal template.
        config = json.loads(journal_path.read_text(encoding="utf-8"))
    else:
        # Load one of the templates shipped with TARDIS.
        journal_name = str(journal).lower()
        if journal_name not in BUILTIN_JOURNALS:
            supported = ", ".join(BUILTIN_JOURNALS)
            msg = (
                f"Unsupported journal {journal!r}; choose one of: "
                f"{supported}, or provide a JSON configuration path"
            )
            raise ValueError(msg)
        templates = json.loads(
            files("tardis.data")
            .joinpath(JOURNAL_TEMPLATES_FILE)
            .read_text(encoding="utf-8")
        )
        config = templates[journal_name]

    # Check the two settings needed to render the table.
    if not isinstance(config.get("bold_total"), bool):
        msg = "Journal configuration 'bold_total' must be a boolean"
        raise ValueError(msg)
    template = config.get("template")
    if not isinstance(template, str) or "{{TABLE_ROWS}}" not in template:
        msg = "Journal configuration template must contain '{{TABLE_ROWS}}'"
        raise ValueError(msg)
    return config


def _render_atom_data_summary(
    summary: pd.DataFrame,
    config: dict[str, object],
) -> str:
    """Render an atomic-data summary with a journal configuration."""
    rows = []

    # Turn each summary record into one LaTeX table row.
    for record in summary.itertuples(index=False):
        if record.Element == "Total" and config["bold_total"]:
            row = (
                rf"\textbf{{Total}} &  & \textbf{{{record.Levels}}} & "
                rf"\textbf{{{record.Lines}}} \\"
            )
        else:
            row = (
                f"{record.Element} & {record[1]} & "
                f"{record.Levels} & {record.Lines} " + r"\\"
            )
        rows.append(row)

    # Place the finished rows inside the selected journal template.
    return str(config["template"]).replace("{{TABLE_ROWS}}", "\n".join(rows))

# The user-facing public function.
# Read tables in, build summary, load format, output .tex and .txt file, return dataframe.
def export_atom_data_summary(
    atom_data: AtomData | str | Path,
    output_path: str | Path,
    *,
    journal: str | Path = "aas",
    elements: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Export a publication-ready atomic-data summary table.

    Parameters
    ----------
    atom_data : tardis.io.atom_data.AtomData or str or pathlib.Path
        Loaded TARDIS atomic data or a Carsus/TARDIS atomic HDF file.
    output_path : str or pathlib.Path
        Output basename or ``.tex`` path. Matching ``.tex`` and ``.txt`` files
        are written.
    journal : {"aa", "aas", "mnras", "nature"} or pathlib.Path, optional
        Built-in journal style or path to a custom JSON configuration.
    elements : collections.abc.Iterable[str], optional
        Element symbols to include. All available elements are included by
        default.

    Returns
    -------
    pandas.DataFrame
        The summary written to disk.
    """
    # Read the atomic tables and build the shared summary data.
    levels, lines = _read_atom_data_tables(atom_data)
    summary = build_atom_data_summary(levels, lines, elements=elements)

    # Load the journal layout that will wrap the LaTeX rows.
    config = _load_journal_config(journal)

    # Treat either a basename or a .tex filename as the output basename.
    output_path = Path(output_path)
    output_base = (
        output_path.with_suffix("") if output_path.suffix else output_path
    )
    # Create the destination folder and write both output formats.
    output_base.parent.mkdir(parents=True, exist_ok=True)
    output_base.with_suffix(".tex").write_text(
        _render_atom_data_summary(summary, config),
        encoding="utf-8",
    )
    output_base.with_suffix(".txt").write_text(
        summary.to_string(index=False) + "\n",
        encoding="utf-8",
    )
    return summary
