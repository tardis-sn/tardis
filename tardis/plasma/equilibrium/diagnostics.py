"""Replay a retained level solution without transport or global plasma solves."""

from collections.abc import Hashable

import numpy as np
import pandas as pd
from astropy import units as u

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.evaluator import (
    PlasmaEquilibriumEvaluator,
    calculate_nlte_level_population_residual,
)


def diagnose_level_closure(
    evaluator: PlasmaEquilibriumEvaluator,
    shell: Hashable,
    *,
    resolve: bool = False,
) -> dict[str, pd.DataFrame]:
    """Compare inner closure with staged changes to the terminal state.

    Parameters
    ----------
    evaluator : PlasmaEquilibriumEvaluator
        Evaluator with a completed ``last_level_snapshot``. Use it before
        changing its fixed scientific inputs. Only its latest completed call
        is retained; a failed subsequent call does not replace that snapshot.
    shell : hashable
        Shell column label, not a positional array index.
    resolve : bool, optional
        Retry only this shell's inner level optimizer using the retained
        fractions. False performs residual checks only. True can be expensive
        for a difficult shell, but never runs the ion solver or transport.

    Returns
    -------
    dict of str to pandas.DataFrame
        ``context`` identifies the call and tentative diagnostic branch;
        ``stages`` reports worst all-row residuals; ``levels`` gives fractions,
        net/gross flows and continuum terms by stage; ``beta`` compares escape
        probabilities; ``rates`` compares trial/final density-specific continuum
        rates. Flows are s^-1 per selected-ion population, not cm^-3 s^-1.

    Notes
    -----
    Stages successively change electron density, absolute opacity geometry,
    represented fractions and the authoritative adjacent-ion ratio. Earlier
    stages eliminate that ratio algebraically, as the inner solver does.
    Differences isolate sensitivity, not independent proof of a code defect.
    The reported terminal residual must agree with its replay before drawing
    conclusions. No population, rate, acceptance or evaluator snapshot is
    modified. Rate preparation uses all shell columns to preserve label/index
    contracts, but residuals and any optional nonlinear retry concern one shell.
    """
    snapshot = getattr(evaluator, "last_level_snapshot", None)
    if snapshot is None:
        raise ValueError(
            "No completed evaluator snapshot is available for diagnosis."
        )
    evaluation = snapshot.evaluation
    shell_idx = evaluation.normalized_population.columns.get_loc(shell)
    species_index = evaluation.normalized_population.index
    fractions = evaluation.normalized_population.iloc[:, shell_idx].to_numpy(
        copy=True
    )
    absolute = evaluation.absolute_level_population.iloc[:, shell_idx].to_numpy(
        copy=True
    )
    represented = evaluation.absolute_level_population.loc[
        species_index, shell
    ].to_numpy(copy=True)
    represented_sum = represented.sum()
    represented_fractions = (
        represented / represented_sum
        if represented_sum > 0
        else fractions.copy()
    )
    trial_densities = evaluation.trial_electron_density
    final_densities = pd.Series(
        np.asarray(
            trial_densities
            if evaluation.charge_solved_electron_density is None
            else evaluation.charge_solved_electron_density
        ),
        index=trial_densities.index,
    )
    trial_density = float(trial_densities.iloc[shell_idx])
    final_density = float(final_densities.iloc[shell_idx])
    temperature = float(snapshot.electron_temperature[shell_idx])
    continuum = snapshot.continuum_rates[shell_idx]
    geometry = snapshot.shell_number_densities[shell_idx]
    sobolev = evaluator.sobolev_inputs[shell_idx]
    bound_trial = snapshot.bound_bound_rates[shell_idx]
    bound_final = bound_trial
    if not np.array_equal(
        trial_densities.to_numpy(), np.asarray(final_densities)
    ):
        bound_final = evaluator.rate_matrix_solver.prepare_species_rates(
            evaluator.j_blues,
            ThermalElectronEnergyDistribution(
                0.0 * u.erg,
                snapshot.electron_temperature * u.K,
                np.asarray(final_densities) / u.cm**3,
            ),
            evaluator.equilibrium_species,
            evaluator.sobolev_inputs[0].line_index,
            None
            if evaluator.fixed_rates is None
            else evaluator.fixed_rates.bound_bound_rates,
        )[shell_idx]
    lower_density = np.nan
    upper_density = np.nan
    authoritative_ratio = None
    if evaluation.ion_population is not None:
        species = evaluator.equilibrium_species
        lower_density = float(evaluation.ion_population.loc[species, shell])
        upper_density = float(
            evaluation.ion_population.loc[(species[0], species[1] + 1), shell]
        )
        if represented_sum > 0.0:
            authoritative_ratio = upper_density / represented_sum

    # Each tuple states precisely which inputs differ; no full evaluate call.
    configurations = [
        ("inner", fractions, trial_density, bound_trial, None, None),
        ("final_density", fractions, final_density, bound_final, None, None),
        (
            "final_geometry",
            fractions,
            final_density,
            bound_final,
            absolute,
            None,
        ),
        (
            "represented_fractions",
            represented_fractions,
            final_density,
            bound_final,
            absolute,
            None,
        ),
        (
            "terminal",
            represented_fractions,
            final_density,
            bound_final,
            absolute,
            authoritative_ratio,
        ),
    ]
    if resolve:
        cache_names = (
            "_level_solution_cache",
            "_level_solution_cache_species",
            "_level_solution_cache_index",
            "_level_solution_cache_columns",
        )
        saved_attributes = {
            name: getattr(evaluator, name)
            for name in cache_names
            if hasattr(evaluator, name)
        }
        cache = saved_attributes.get("_level_solution_cache")
        saved_entries = None if cache is None else cache.copy()
        saved_values = (
            None
            if cache is None
            else {key: value.copy() for key, value in cache.items()}
        )
        try:
            resolved, _, _, _ = evaluator._calculate_level_solution(
                shell_idx,
                trial_density,
                temperature,
                fractions.copy(),
                continuum,
                geometry,
                sobolev,
                bound_trial,
            )
        finally:
            if cache is not None:
                for key, value in saved_entries.items():
                    value[:] = saved_values[key]
                cache.clear()
                cache.update(saved_entries)
            for name in cache_names:
                if name in saved_attributes:
                    setattr(evaluator, name, saved_attributes[name])
                elif hasattr(evaluator, name):
                    delattr(evaluator, name)
        configurations.append(
            ("resolved_inner", resolved, trial_density, bound_trial, None, None)
        )
        resolved_absolute = absolute.copy()
        resolved_absolute[geometry.species_level_positions] = (
            lower_density * resolved
            if np.isfinite(lower_density)
            else geometry.species_number_density
            / (1.0 + evaluation.diagnostic_ion_ratio.iloc[shell_idx])
            * resolved
        )
        resolved_sum = resolved_absolute[geometry.species_level_positions].sum()
        resolved_ratio = (
            upper_density / resolved_sum
            if resolved_sum > 0 and np.isfinite(upper_density)
            else None
        )
        configurations.append(
            (
                "resolved_terminal",
                resolved,
                final_density,
                bound_final,
                resolved_absolute,
                resolved_ratio,
            )
        )

    summaries = []
    level_frames = {}
    beta_frames = {}
    rate_frames = {}
    for name, candidate, density, bound, level_density, ratio in configurations:
        rates = evaluator._calculate_level_equation_rates(
            shell_idx, density, continuum
        )
        recombination = rates.recombination.copy()
        if evaluator.equilibrium_species in evaluator.case_b_species:
            recombination[0] = 0.0
        eliminated_ratio = float(
            rates.ionization @ candidate / recombination.sum()
        )
        used_ratio = eliminated_ratio if ratio is None else ratio
        budget = {}
        residual, beta, _ = calculate_nlte_level_population_residual(
            candidate,
            rates,
            evaluator.rate_matrix_solver,
            None,
            None,
            evaluator.equilibrium_species,
            geometry,
            sobolev,
            level_density=level_density,
            bound_bound_rates=bound,
            case_b=evaluator.equilibrium_species in evaluator.case_b_species,
            normalize_residual=True,
            relative_tolerance=evaluator.level_residual_tolerance,
            continuum_population_ratio=ratio,
            replace_normalization_row=False,
            diagnostic_budget=budget,
        )
        beta_frames[name] = beta
        rate_frames[name] = pd.DataFrame(
            {"ionization": rates.ionization, "recombination": recombination},
            index=species_index,
        )
        level_frames[name] = pd.DataFrame(
            {
                "fraction": candidate,
                "relative_residual": residual,
                "net_flow": budget["net_flow"],
                "gross_flow": budget["gross_flow"],
                "flow_scale": budget["flow_scale"],
                "ionization_loss": rates.ionization * candidate,
                "recombination_source": recombination * used_ratio,
                "bound_bound_net_flow": (
                    budget["net_flow"]
                    + rates.ionization * candidate
                    - recombination * used_ratio
                ),
            },
            index=species_index,
        )
        worst_idx = int(np.argmax(np.abs(residual)))
        summaries.append(
            {
                "stage": name,
                "electron_density": density,
                "eliminated_ion_ratio": eliminated_ratio,
                "used_ion_ratio": used_ratio,
                "total_ionization": float(rates.ionization @ candidate),
                "total_recombination": float(recombination.sum() * used_ratio),
                "max_relative_residual": float(np.max(np.abs(residual))),
                "worst_level": species_index[worst_idx],
                "max_fraction_change_from_inner": float(
                    np.max(np.abs(candidate - fractions))
                ),
                "max_beta_change_from_inner": float(
                    np.max(np.abs(beta - beta_frames["inner"]), initial=0.0)
                ),
            }
        )
    stages = pd.DataFrame(summaries).set_index("stage")
    reported = evaluation.relative_level_residual
    agrees = reported is not None and np.allclose(
        level_frames["terminal"]["relative_residual"],
        reported[shell],
        rtol=64 * np.finfo(float).eps,
        atol=64 * np.finfo(float).eps,
    )
    inner_closed = (
        stages.loc["inner", "max_relative_residual"]
        <= evaluator.level_residual_tolerance
    )
    terminal_closed = (
        stages.loc["terminal", "max_relative_residual"]
        <= evaluator.level_residual_tolerance
    )
    branch = "closed" if terminal_closed else "post_level_state_changed"
    if not inner_closed:
        branch = "inner_level_solve_open"
    if not agrees:
        branch = "replay_context_mismatch"
    elif represented_sum == 0.0 and evaluation.ion_population is not None:
        # The legacy terminal calculation falls back to an eliminated ratio
        # here. Do not interpret that as verified authoritative ion closure.
        branch = "zero_selected_density"
    context = pd.DataFrame(
        {
            "evaluation_count": snapshot.evaluation_count,
            "species": [evaluator.equilibrium_species],
            "electron_temperature": temperature,
            "selected_ion_density": lower_density,
            "upper_ion_density": upper_density,
            "represented_level_sum": represented_sum,
            "inner_pair_density": geometry.species_number_density,
            "relative_tolerance": evaluator.level_residual_tolerance,
            "terminal_replay_matches_reported": bool(agrees),
            "authoritative_ratio_available": authoritative_ratio is not None,
            "absolute_insignificance_accepted": (
                False
                if evaluation.absolute_level_budget is None
                else bool(
                    evaluation.absolute_level_budget.loc[shell, "accepted"]
                )
            ),
            "branch": branch,
        },
        index=pd.Index([shell], name="shell"),
    )
    return {
        "context": context,
        "stages": stages,
        "levels": pd.concat(level_frames, axis=1, names=["stage", "quantity"]),
        "beta": pd.DataFrame(
            {
                **beta_frames,
                "reported_inner": evaluation.trial_beta_sobolev.iloc[
                    :, shell_idx
                ].to_numpy(),
                "reported_terminal": evaluation.beta_sobolev.iloc[
                    :, shell_idx
                ].to_numpy(),
            },
            index=sobolev.line_index,
        ),
        "rates": pd.concat(rate_frames, axis=1, names=["stage", "quantity"]),
    }
