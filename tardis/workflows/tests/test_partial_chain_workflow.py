"""Comprehensive two-iteration verification of Type IIP reduction modes."""

from pathlib import Path

import numpy as np
import pytest

from tardis.io.configuration.config_reader import Configuration
from tardis.plasma.equilibrium.thermal_balance import (
    ElectronThermalBalanceSolver,
)
from tardis.workflows.type_iip_workflow import TypeIIPWorkflow


@pytest.fixture
def partial_chain_config(tardis_regression_path: Path) -> Configuration:
    """Use the existing mixed H/He network, without changing enabled physics."""
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/iip_hydrogen_helium.yml"
    )
    config.atom_data = (
        tardis_regression_path
        / "atom_data/christians_atomdata_converted_04Dec25.h5"
    )
    config.montecarlo.nthreads = 1
    return config


@pytest.mark.parametrize(
    "mode,elements,include_thermal_pool,include_ionization_pool",
    [
        ("full", [], False, False),
        ("partial", ["H"], False, False),
        ("partial", ["H"], True, False),
        ("partial", ["H"], False, True),
        ("partial", ["H"], True, True),
        ("partial", ["He"], False, False),
        ("partial", ["H", "He"], True, True),
        ("explicit", [], False, False),
    ],
)
def test_coupled_workflow_thermal_closure(
    partial_chain_config: Configuration,
    mode: str,
    elements: list[str],
    include_thermal_pool: bool,
    include_ionization_pool: bool,
) -> None:
    """Each mode closes the same thermal equations through two transport steps.

    Regime: the fixture's 200000 packets, mixed H/He, one plasma update and
    final transport. Retain the physical fixture's packet count because smaller
    development samples can trigger the known nonlinear thermal-solver basin.
    Independent verification: the accepted physical balance residual is tested
    against the solver's predeclared closure criteria, rather than just finite
    output. Cross-mode stochastic comparisons are separate validation evidence.
    """
    config = partial_chain_config
    config.plasma.continuum_interaction.markov_chain.mode = mode
    config.plasma.continuum_interaction.markov_chain.elements = elements
    config.plasma.continuum_interaction.markov_chain.include_thermal_pool = (
        include_thermal_pool
    )
    config.plasma.continuum_interaction.markov_chain.include_ionization_pool = (
        include_ionization_pool
    )
    workflow = TypeIIPWorkflow(config)
    workflow.show_progress_bars = False
    workflow.run()
    ElectronThermalBalanceSolver.validate_evaluation(
        workflow._thermal_balance_evaluation
    )
    assert workflow.completed_iterations == 1
    solver = workflow.macro_atom_solver
    metadata = solver.computed_metadata[0]
    continuum = workflow.transport_state.continuum_state_numba
    np.testing.assert_array_equal(
        continuum.macro_level_to_reduced_idx, solver.macro_level_to_reduced_idx
    )
    np.testing.assert_array_equal(
        continuum.reduced_to_macro_level_idx, solver.reduced_to_macro_level_idx
    )
    assert len(continuum.photo_ion_activation_idx) == len(
        workflow.continuum_opacity_state.level2continuum_idx
    )
    photo_ion_block = metadata.loc[
        metadata.source == ("i", -99, -99), "source_level_idx"
    ].unique()[0]
    assert np.all(continuum.photo_ion_activation_idx == photo_ion_block)
    for name in (
        "photo_ion_estimator",
        "photo_ion_estimator_statistics",
        "stim_recomb_estimator",
        "bf_heating_estimator",
        "stim_recomb_cooling_estimator",
    ):
        estimator = getattr(workflow.transport_state.estimators_continuum, name)
        assert estimator.shape == (
            len(workflow.continuum_opacity_state.level2continuum_idx),
            workflow.simulation_state.geometry.no_of_shells_active,
        )
    if mode == "partial":
        selected = set(continuum.reduced_to_macro_level_idx)
        for pool, included in (
            ("k", include_thermal_pool),
            ("i", include_ionization_pool),
        ):
            pool_levels = set(
                metadata.loc[
                    metadata.source.map(lambda level: level[0] == pool),
                    "source_level_idx",
                ]
            )
            if included:
                assert pool_levels <= selected
            else:
                assert pool_levels.isdisjoint(selected)
