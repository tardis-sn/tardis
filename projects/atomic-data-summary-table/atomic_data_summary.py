"""Project copy of atomic_data_summary for PR packaging.

This is a repo-local copy placed under projects/atomic-data-summary-table
so it can be included in the PR without modifying package layout.
"""

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

def _format_ion_stage_range(ion_numbers: Iterable[int]) -> str:
    stages = sorted({int(ion_number) + 1 for ion_number in ion_numbers})
    runs = []
    start = previous = stages[0]
    for stage in stages[1:]:
        if stage != previous + 1:
            runs.append((start, previous))
            start = stage
        previous = stage
    runs.append((start, previous))
    return ", ".join(
        int_to_roman(start)
        if start == end
        else f"{int_to_roman(start)}--{int_to_roman(end)}"
        for start, end in runs
    )


def _get_index_columns(dataframe: pd.DataFrame) -> pd.DataFrame:
    required_columns = {"atomic_number", "ion_number"}
    if required_columns.issubset(dataframe.columns):
        return dataframe
    return dataframe.reset_index()


def build_atom_data_summary(
    levels: pd.DataFrame,
    lines: pd.DataFrame,
    elements: Iterable[str] | None = None,
) -> pd.DataFrame:
    levels = _get_index_columns(levels)
    lines = _get_index_columns(lines)
    grouping_columns = ["atomic_number", "ion_number"]
    counts_by_ion = pd.concat(
        [
            levels.groupby(grouping_columns).size().rename("Levels"),
            lines.groupby(grouping_columns).size().rename("Lines"),
        ],
        axis=1,
    )
    eligible_ions = counts_by_ion[(counts_by_ion["Levels"] > 1) & counts_by_ion["Lines"].notna()].reset_index()
    if elements is not None:
        atomic_numbers = {element_symbol2atomic_number(element) for element in elements}
        eligible_ions = eligible_ions[eligible_ions["atomic_number"].isin(atomic_numbers)]
    grouped_ions = eligible_ions.groupby("atomic_number", sort=True)
    summary = pd.concat(
        [
            grouped_ions["ion_number"].apply(_format_ion_stage_range).rename("Ion stages"),
            grouped_ions[["Levels", "Lines"]].sum(),
        ],
        axis=1,
    )
    summary.insert(
        0,
        "Element",
        [
            rf"$\mathrm{{{atomic_number2element_symbol(atomic_number)}}}"
            rf"_{{{atomic_number}}}$"
            for atomic_number in summary.index
        ],
    )
    summary = summary.reset_index(drop=True)
    summary = summary.rename(columns={"Levels": "No. of Levels", "Lines": "No. of Lines"})
    summary[["No. of Levels", "No. of Lines"]] = summary[["No. of Levels", "No. of Lines"]].astype(int)
    summary.loc[len(summary)] = [
        "Total",
        "",
        int(summary["No. of Levels"].sum()),
        int(summary["No. of Lines"].sum()),
    ]
    return summary


def _read_atom_data_tables(atom_data: AtomData | str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    if isinstance(atom_data, (str, Path)):
        with pd.HDFStore(atom_data, mode="r") as store:
            return store["levels_data"], store["lines_data"]
    return atom_data.levels, atom_data.lines


def _load_journal_config(journal: str | Path) -> dict[str, object]:
    journal_path = Path(journal)
    if journal_path.suffix.lower() == ".json":
        config = json.loads(journal_path.read_text(encoding="utf-8"))
    else:
        journal_name = str(journal).lower()
        if journal_name not in BUILTIN_JOURNALS:
            supported = ", ".join(BUILTIN_JOURNALS)
            raise ValueError(f"Unsupported journal {journal!r}; choose one of: {supported}")
        templates = json.loads(files("tardis.data").joinpath(JOURNAL_TEMPLATES_FILE).read_text(encoding="utf-8"))
        config = templates[journal_name]
    if not isinstance(config.get("bold_total"), bool):
        raise ValueError("Journal configuration 'bold_total' must be a boolean")
    template = config.get("template")
    if not isinstance(template, str) or "{{TABLE_ROWS}}" not in template:
        raise ValueError("Journal configuration template must contain '{{TABLE_ROWS}}'")
    return config


def _render_atom_data_summary(summary: pd.DataFrame, config: dict[str, object]) -> str:
    rows = []
    for record in summary.itertuples(index=False):
        elem = record[0]
        ion_stage = record[1]
        n_levels = record[2]
        n_lines = record[3]
        if elem == "Total":
            rows.append(r"\\hline")
            if config["bold_total"]:
                row = (rf"\\textbf{{Total}} &  & \\textbf{{{n_levels}}} & " rf"\\textbf{{{n_lines}}} \\\")
            else:
                row = f"{elem} & {ion_stage} & {n_levels} & {n_lines} " + r"\\"
        else:
            row = f"{elem} & {ion_stage} & {n_levels} & {n_lines} " + r"\\"
        rows.append(row)
    return str(config["template"]).replace("{{TABLE_ROWS}}", "\n".join(rows))


def export_atom_data_summary(atom_data: AtomData | str | Path, output_path: str | Path, journal: str | Path = "aas", elements: Iterable[str] | None = None) -> pd.DataFrame:
    levels, lines = _read_atom_data_tables(atom_data)
    summary = build_atom_data_summary(levels, lines, elements=elements)
    config = _load_journal_config(journal)
    output_path = Path(output_path)
    output_base = (output_path.with_suffix("") if output_path.suffix else output_path)
    output_base.parent.mkdir(parents=True, exist_ok=True)
    output_base.with_suffix(".tex").write_text(_render_atom_data_summary(summary, config), encoding="utf-8")
    output_base.with_suffix(".txt").write_text(summary.to_string(index=False) + "\n", encoding="utf-8")
    return summary
