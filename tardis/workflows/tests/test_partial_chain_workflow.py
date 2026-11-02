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
    "mode,elements",
    [
        ("full", []),
        ("partial", ["H"]),
        ("partial", ["He"]),
        ("partial", ["H", "He"]),
        ("explicit", []),
    ],
)
def test_coupled_workflow_thermal_closure(
    partial_chain_config: Configuration, mode: str, elements: list[str]
) -> None:
    """Each mode closes the same thermal equations through two transport steps.

    Regime: 20000 packets, mixed H/He, one plasma update and final transport.
    Independent verification: the accepted physical balance residual is tested
    against the solver's predeclared closure criteria, rather than just finite
    output. Cross-mode stochastic comparisons are separate validation evidence.
    """
    config = partial_chain_config
    config.montecarlo.no_of_packets = 20000
    config.montecarlo.last_no_of_packets = 20000
    config.plasma.continuum_interaction.markov_chain.mode = mode
    config.plasma.continuum_interaction.markov_chain.elements = elements
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
        pool_levels = metadata.loc[
            metadata.source.map(lambda level: level[0] in ("i", "k")),
            "source_level_idx",
        ]
        assert not set(continuum.reduced_to_macro_level_idx).intersection(
            pool_levels
        )
