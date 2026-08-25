from types import SimpleNamespace

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal, assert_series_equal
from scipy.integrate import trapezoid

from tardis import constants as const
from tardis.iip_plasma.continuum.base import ContinuumProcess
from tardis.iip_plasma.properties.continuum import (
    IIpWorkflowContinuumConnectors,
    integrate_array_by_level_groups,
)


def test_normalize_transition_probabilities_matches_groupby_transform():
    """Compare transition normalization with the pandas groupby baseline."""
    index = pd.MultiIndex.from_tuples(
        [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5)],
        names=["source_level_idx", "destination_level_idx"],
    )
    transition_probabilities = pd.DataFrame(
        {
            "transition_type": [1, 0, 2, -1, 0],
            "lines_idx": [10, 11, -1, 12, 13],
            0: [1.0, 3.0, 0.0, 0.0, 5.0],
            1: [2.0, 2.0, 0.0, 0.0, 15.0],
        },
        index=index,
    )
    expected_probabilities = (
        transition_probabilities.iloc[:, 2:]
        .groupby(level=0)
        .transform(lambda values: values / values.sum())
    )
    expected = pd.concat(
        [transition_probabilities.iloc[:, :2], expected_probabilities],
        axis=1,
    ).fillna(0.0)

    actual = ContinuumProcess._normalize_transition_probabilities(
        transition_probabilities, no_ref_columns=2
    )

    assert_frame_equal(actual, expected)
    assert np.issubdtype(actual["transition_type"].dtype, np.integer)
    assert np.issubdtype(actual["lines_idx"].dtype, np.integer)


def test_integrate_array_by_level_groups_matches_groupby_apply():
    index = pd.MultiIndex.from_tuples(
        [
            (1, 0, 0),
            (1, 0, 0),
            (1, 0, 0),
            (2, 1, 0),
            (2, 1, 0),
            (2, 1, 0),
            (2, 1, 1),
            (2, 1, 1),
        ],
        names=["atomic_number", "ion_number", "level_number"],
    )
    nu = pd.Series([1.0, 1.5, 2.0, 3.0, 3.25, 4.0, 5.0, 6.0], index=index)
    values = np.array(
        [
            [1.0, 2.0, 3.0],
            [1.5, 2.5, 3.5],
            [2.0, 3.0, 4.0],
            [0.5, 1.0, 1.5],
            [0.75, 1.25, 1.75],
            [1.0, 1.5, 2.0],
            [1.25, 1.75, 2.25],
            [1.5, 2.0, 2.5],
        ]
    )

    actual = integrate_array_by_level_groups(values, nu)

    frame = pd.DataFrame(values, index=index)
    frame.insert(0, "nu", nu)
    grouped = frame.groupby(level=[0, 1, 2])
    expected = pd.DataFrame(
        {
            shell: grouped.apply(
                lambda sub, shell=shell: trapezoid(sub[shell], sub["nu"])
            )
            for shell in range(values.shape[1])
        }
    )

    assert_frame_equal(actual, expected)


def test_integrate_array_by_level_groups_matches_series_groupby_apply():
    index = pd.MultiIndex.from_tuples(
        [
            (1, 0, 0),
            (1, 0, 0),
            (1, 0, 0),
            (2, 1, 0),
            (2, 1, 0),
        ],
        names=["atomic_number", "ion_number", "level_number"],
    )
    nu = pd.Series([1.0, 1.5, 2.0, 3.0, 4.0], index=index)
    values = np.array([1.0, 1.5, 2.0, 0.5, 1.0])

    level_groups = tuple(nu.groupby(level=[0, 1, 2]).indices.items())
    actual = integrate_array_by_level_groups(values, nu, level_groups)

    frame = pd.DataFrame({"values": values, "nu": nu}, index=index)
    expected = frame.groupby(level=[0, 1, 2]).apply(
        lambda sub: trapezoid(sub["values"], sub["nu"])
    )

    assert_series_equal(actual, expected)


def test_iip_bound_free_opacity_remains_finite_after_lte_population_cutoff() -> None:
    """Calculate finite opacity after outer-shell LTE populations are clipped."""
    level_names = ["atomic_number", "ion_number", "level_number"]
    level_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1), (1, 1, 0)], names=level_names
    )
    photoionization_index = pd.MultiIndex.from_tuples(
        [(1, 0, 1), (1, 0, 1)], names=level_names
    )
    photoionization_data = pd.DataFrame(
        {
            "nu": [1.0e15, 2.0e15],
            "x_sect": [1.0e-18, 2.0e-18],
        },
        index=photoionization_index,
    )
    atomic_data = SimpleNamespace(
        photoionization_data=photoionization_data,
        macro_atom_references=pd.DataFrame(
            {"references_idx": [0, 1, 2]}, index=level_index
        ),
        levels=pd.DataFrame({"energy": [0.0, 1.0, 2.0]}, index=level_index),
    )
    continuum_index = pd.MultiIndex.from_tuples(
        [(1, 0, 1)], names=level_names
    )
    ion_index = pd.MultiIndex.from_tuples(
        [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
    )
    transition_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0, 1)],
        names=[
            "atomic_number",
            "ion_number",
            "level_number_lower",
            "level_number_upper",
        ],
    )
    electron_densities = np.ones(2)
    ion_number_density = pd.DataFrame(
        [[1.0, 1.0], [1.0, 0.0]], index=ion_index
    )
    t_electrons = np.full(2, 1.0e10)
    level_number_density = pd.DataFrame(
        [[1.0, 1.0], [2.0, 0.5], [1.0, 1.0]], index=level_index
    )
    phi_lucy = pd.DataFrame([[0.5, 0.25]], index=continuum_index)
    plasma_inputs = {
        "atomic_data": atomic_data,
        "alpha_sp": pd.DataFrame([[1.0, 1.0]], index=continuum_index),
        "electron_densities": electron_densities,
        "ion_number_density": ion_number_density,
        "lte_ion_number_density": pd.DataFrame(
            [[1.0, 0.0], [1.0, 0.0]], index=ion_index
        ),
        "t_electrons": t_electrons,
        "level_number_density": level_number_density,
        "lte_level_number_density": pd.DataFrame(
            [[1.0, 0.0], [0.5, 0.0], [1.0, 0.0]], index=level_index
        ),
        "phi_lucy": phi_lucy,
        "gamma": None,
        "coll_deexc_coeff": pd.DataFrame(
            [[1.0, 1.0]], index=transition_index
        ),
    }
    continuum_connectors = IIpWorkflowContinuumConnectors(
        SimpleNamespace(get_value=plasma_inputs.__getitem__)
    )

    continuum_connectors.update()

    boltzmann_factor = np.exp(
        -photoionization_data["nu"].values[:, np.newaxis]
        / t_electrons
        * (const.h.cgs.value / const.k_B.cgs.value)
    )
    expected = (
        level_number_density.loc[photoionization_index].values
        - phi_lucy.loc[photoionization_index].values
        * ion_number_density.loc[[(1, 1), (1, 1)]].values
        * electron_densities
        * boltzmann_factor
    )
    expected *= photoionization_data["x_sect"].values[:, np.newaxis]

    np.testing.assert_allclose(continuum_connectors.chi_bf.values, expected)
    assert np.isfinite(continuum_connectors.chi_bf.values).all()
    assert np.isfinite(continuum_connectors.p_fb_deactivation.values).all()
    np.testing.assert_array_equal(
        continuum_connectors.p_fb_deactivation.iloc[:, 1].values,
        np.zeros(1),
    )
