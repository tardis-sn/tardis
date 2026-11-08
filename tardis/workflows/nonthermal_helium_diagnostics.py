"""Inspect retained helium coupling states without running a plasma solver.

This module imports no TARDIS modules, so it can also be loaded with
``runpy.run_path`` in a notebook that still holds an older workflow instance.
"""

import numpy as np
import pandas as pd


def canonicalize_underflowed_populations(
    levels: pd.DataFrame, ions: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Zero orphan subnormal parents whose represented levels are exactly zero.

    Parameters
    ----------
    levels : pandas.DataFrame
        Level number densities in cm^-3, with atomic/ion/level row labels.
    ions : pandas.DataFrame
        Parent ion densities in cm^-3, with atomic/ion row labels.

    Returns
    -------
    levels, ions : tuple of pandas.DataFrame
        Populations with only allocation-underflowed parents zeroed, without
        changing caller-owned tables. The level
        table is unchanged, since every affected level is already zero. The
        ion table is copied only when a correction is necessary.

    Notes
    -----
    Binary64 subnormals lie below `finfo(float64).tiny`, the smallest normal
    number, and extend down to `nextafter(0, 1)`. Under round-to-nearest-even,
    a level share rounds to zero only up to half that smallest subnormal.
    Thus normalized nonnegative allocation among N represented levels can
    vanish entirely only for a parent at most floor(N/2) subnormal quanta.
    Only an all-zero distribution within that bound is canonicalized; nonzero
    subnormal levels, missing parents and contaminated data are not repaired.
    This is representability handling, not an astrophysical abundance floor.
    """
    # Nullable missing values are not evidence of an exactly zero level.
    level_groups = (
        levels.eq(0.0)
        .fillna(False)
        .groupby(level=["atomic_number", "ion_number"])
    )
    zero_levels = level_groups.all()
    # Do not halve the smallest subnormal first: that already rounds to zero.
    allocation_bound = (level_groups.size() // 2) * np.nextafter(0.0, 1.0)
    aligned_ions = ions.reindex(
        index=zero_levels.index, columns=zero_levels.columns
    )
    underflowed = (
        (
            zero_levels
            & (aligned_ions > 0.0)
            & (aligned_ions < np.finfo(np.float64).tiny)
            & aligned_ions.le(allocation_bound, axis="index")
        )
        .reindex(index=ions.index, columns=ions.columns, fill_value=False)
        .fillna(False)
    )
    if underflowed.to_numpy().any():
        ions = ions.mask(underflowed, 0.0)
    return levels, ions


def diagnose_level_relaxation(
    current_levels: pd.DataFrame,
    candidate_levels: pd.DataFrame,
    current_ions: pd.DataFrame,
    candidate_ions: pd.DataFrame,
    relaxation: float = 0.5,
    ion_relaxation_factors: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Audit a selective level update using snapshots, without changing them.

    Parameters
    ----------
    current_levels, candidate_levels : pandas.DataFrame
        Input and unrelaxed candidate level densities in cm^-3, indexed by
        atomic number, ion number, and level number; columns are shell labels.
    current_ions, candidate_ions : pandas.DataFrame
        Corresponding ion densities in cm^-3, indexed by atomic and ion number.
    relaxation : float, optional
        Fixed ion-population damping fraction when explicit factors are absent.
    ion_relaxation_factors : pandas.DataFrame, optional
        Actual per-element factors, indexed by atomic number with shell
        columns. Supply these when auditing an adaptive coupling update.

    Returns
    -------
    pandas.DataFrame
        One row per represented candidate ion and shell. Reports parent-ion
        availability, level sums and nonfinite counts in each snapshot, the
        predicted nonfinite count after selective relaxation, and whether that
        update creates nonfinite levels from finite inputs. ``classification``
        distinguishes missing parents, contaminated inputs/candidates, and
        finite zero distributions with a positive ion total. Input/candidate
        columns describe the raw snapshots; ``*_underflow_repaired`` flags
        parents zeroed for the projected production update. Sums ending in
        ``skipna`` deliberately reproduce pandas' all-missing-to-zero sum;
        ``strict`` sums remain NaN if any level is nonfinite. This audit can
        identify a propagation mechanism, not reconstruct the first corrupt
        iteration from an already contaminated snapshot.
    """
    ion_names = ["atomic_number", "ion_number"]
    columns = candidate_levels.columns
    current_levels = current_levels.reindex(
        index=candidate_levels.index, columns=columns
    )
    current_ions = current_ions.reindex(columns=columns)
    candidate_ions = candidate_ions.reindex(columns=columns)
    groups = candidate_levels.groupby(level=ion_names).sum().index
    level_groups = candidate_levels.index.droplevel("level_number")
    report_index = pd.MultiIndex.from_tuples(
        [(*ion, shell) for ion in groups for shell in columns],
        names=[*ion_names, "shell"],
    )
    report = pd.DataFrame(index=report_index)
    totals = {}
    for name, levels, ions in (
        ("input", current_levels, current_ions),
        ("candidate", candidate_levels, candidate_ions),
    ):
        finite = pd.DataFrame(
            np.isfinite(levels.to_numpy(dtype=float)),
            index=levels.index,
            columns=columns,
        )
        nonfinite_count = (~finite).groupby(level=ion_names).sum()
        level_sums = levels.groupby(level=ion_names).sum()
        totals[name] = levels.groupby(level=ion_names).transform("sum")
        aligned_ions = ions.reindex(index=groups, columns=columns)
        report[f"{name}_nonfinite_level_count"] = (
            nonfinite_count.to_numpy().ravel()
        )
        report[f"{name}_level_sum_skipna"] = level_sums.to_numpy().ravel()
        report[f"{name}_level_sum_strict"] = (
            level_sums.mask(nonfinite_count > 0).to_numpy().ravel()
        )
        report[f"{name}_ion_density"] = aligned_ions.to_numpy().ravel()
        report[f"{name}_parent_present"] = np.repeat(
            groups.isin(ions.index), len(columns)
        )

    current_levels, current_ions = canonicalize_underflowed_populations(
        current_levels, current_ions
    )
    candidate_levels, candidate_ions = canonicalize_underflowed_populations(
        candidate_levels, candidate_ions
    )
    for name, levels, ions in (
        ("input", current_levels, current_ions),
        ("candidate", candidate_levels, candidate_ions),
    ):
        sanitized_ions = ions.reindex(index=groups, columns=columns)
        report[f"{name}_underflow_repaired"] = (
            report[f"{name}_ion_density"] > 0.0
        ) & (sanitized_ions.to_numpy().ravel() == 0.0)
        totals[name] = levels.groupby(level=ion_names).transform("sum")
    relaxed_ions = (
        relaxation * candidate_ions + (1.0 - relaxation) * current_ions
    )
    if ion_relaxation_factors is not None:
        ion_index = relaxed_ions.index
        factors = pd.DataFrame(
            ion_relaxation_factors.reindex(
                index=ion_index.get_level_values("atomic_number"),
                columns=columns,
            ).to_numpy(),
            index=ion_index,
            columns=columns,
        )
        relaxed_ions = factors * candidate_ions + (1 - factors) * current_ions
    fractions = candidate_levels.div(
        totals["candidate"].where(totals["candidate"] != 0.0)
    ).mask(
        totals["candidate"] == 0.0,
        current_levels.div(totals["input"].where(totals["input"] != 0.0)),
    )
    aligned_relaxed_ions = pd.DataFrame(
        relaxed_ions.reindex(index=level_groups, columns=columns).to_numpy(),
        index=candidate_levels.index,
        columns=columns,
    )
    projected = (fractions * aligned_relaxed_ions).mask(
        aligned_relaxed_ions == 0.0, 0.0
    )
    projected, canonical_relaxed_ions = canonicalize_underflowed_populations(
        projected, relaxed_ions
    )
    report["relaxed_underflow_repaired"] = (
        (
            (relaxed_ions.reindex(index=groups, columns=columns) > 0.0)
            & (
                canonical_relaxed_ions.reindex(index=groups, columns=columns)
                == 0.0
            )
        )
        .to_numpy()
        .ravel()
    )
    report["relaxed_ion_density"] = (
        canonical_relaxed_ions.reindex(index=groups, columns=columns)
        .to_numpy()
        .ravel()
    )
    projected_bad = pd.DataFrame(
        ~np.isfinite(projected.to_numpy(dtype=float)),
        index=candidate_levels.index,
        columns=columns,
    )
    newly_bad = (
        projected_bad
        & np.isfinite(current_levels)
        & np.isfinite(candidate_levels)
    )
    report["projected_nonfinite_level_count"] = (
        projected_bad.groupby(level=ion_names).sum().to_numpy().ravel()
    )
    report["introduced_nonfinite_level_count"] = (
        newly_bad.groupby(level=ion_names).sum().to_numpy().ravel()
    )
    report["classification"] = np.select(
        [
            ~report["input_parent_present"]
            | ~report["candidate_parent_present"],
            ~np.isfinite(report["input_ion_density"])
            | ~np.isfinite(report["candidate_ion_density"]),
            (report["input_nonfinite_level_count"] > 0)
            & (report["candidate_nonfinite_level_count"] > 0),
            report["candidate_nonfinite_level_count"] > 0,
            report["input_nonfinite_level_count"] > 0,
            (report["input_level_sum_skipna"] == 0.0)
            & (report["candidate_level_sum_skipna"] == 0.0)
            & (report["relaxed_ion_density"] > 0.0),
            report["input_underflow_repaired"]
            | report["candidate_underflow_repaired"]
            | report["relaxed_underflow_repaired"],
        ],
        [
            "missing_parent_ion",
            "nonfinite_parent_density",
            "nonfinite_input_and_candidate",
            "nonfinite_candidate",
            "nonfinite_input",
            "positive_ion_without_level_distribution",
            "underflow_repaired",
        ],
        default="finite",
    )
    report["affected"] = (report["classification"] != "finite") | (
        report["projected_nonfinite_level_count"] > 0
    )
    return report
