import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis import constants as const
from tardis.plasma.base import BasePlasma
from tardis.plasma.exceptions import IncompleteAtomicData
from tardis.plasma.properties.atomic import IonizationData, Levels, Lines
from tardis.plasma.properties.partition_function import PartitionFunction
from tardis.plasma.properties.radiative_properties import StimulatedEmissionFactor


ionization_data_property = IonizationData(plasma_parent=BasePlasma)
lines_property = Lines(plasma_parent=BasePlasma)
levels_property = Levels(plasma_parent=BasePlasma)
partition_function_property = PartitionFunction(plasma_parent=BasePlasma)




def test_ionization_data_calculate_atomic_property(
    atomic_dataset, selected_atoms, regression_data
):
    actual_ionization_data = ionization_data_property.calculate(
        atomic_dataset, selected_atoms
    )
    expected_ionization_data = regression_data.sync_dataframe(
        actual_ionization_data, key="ionization_data"
    )
    pdt.assert_series_equal(
        actual_ionization_data, expected_ionization_data, atol=0, rtol=1e-15
    )


@pytest.mark.xfail(raises=IncompleteAtomicData)
def test_ionization_data_incomplete_atomic_data(selected_atoms):
    index = pd.MultiIndex.from_tuples(
        [(1, 1), (2, 1)],
        names=["atomic_number", "ion_number"],
    )
    ionization_data = pd.Series(
        [1.0, 2.0], index=index, name="ionization_energy"
    )

    ionization_data_property._filter_atomic_property(
        ionization_data, selected_atoms
    )

def test_levels_calculate(atomic_dataset, selected_atoms, regression_data):
    levels_index, excitation_energy, metastability, g = levels_property.calculate(
        atomic_dataset, selected_atoms
    )

    actual_levels = pd.DataFrame(
        {
            "excitation_energy": excitation_energy,
            "metastability": metastability,
            "g": g,
        },
        index=levels_index,
    )
    expected_levels = regression_data.sync_dataframe(
        actual_levels, key="calculated_levels"
    )
    pdt.assert_frame_equal(
        actual_levels, expected_levels, atol=0, rtol=1e-15
    )


def test_lines_calculate(atomic_dataset, selected_atoms, regression_data):
    actual_lines, _, _, _ = lines_property.calculate(
        atomic_dataset, selected_atoms
    )
    expected_lines = regression_data.sync_dataframe(
        actual_lines, key="calculated_lines"
    )
    pdt.assert_frame_equal(
        actual_lines, expected_lines, atol=0, rtol=1e-15
    )

def test_partition_function_calculate(level_boltzmann_factor_lte, regression_data):
    actual_partition_function_calculate = partition_function_property.calculate(
        level_boltzmann_factor_lte
    )
    expected_partition_function_calculate = regression_data.sync_dataframe(
        actual_partition_function_calculate, key="calculated_partition_function"
    )
    pdt.assert_frame_equal(
        actual_partition_function_calculate, expected_partition_function_calculate, atol=0, rtol=1e-15
    )



@pytest.mark.parametrize(
    "boltzmann_factors, expected_sum",
    [
        ([1.0, 0.5], 1.5),        # two levels — typical case
        ([2.0, 1.0, 0.5], 3.5),   # three levels
        ([1.0], 1.0),             # single level — ground state only
    ],
)
def test_partition_function_sums_boltzmann_factors_within_ion(
    boltzmann_factors, expected_sum
):
    index = pd.MultiIndex.from_tuples(
        [(1, 0, k) for k in range(len(boltzmann_factors))],
        names=["atomic_number", "ion_number", "level_number"],
    )
    level_boltzmann_factor = pd.DataFrame(boltzmann_factors, index=index, columns=[0])
    result = PartitionFunction(None).calculate(level_boltzmann_factor)
    np.testing.assert_allclose(result.loc[(1, 0), 0], expected_sum)


# ── StimulatedEmissionFactor ──────────────────────────────────────────────────


def test_stimulated_emission_factor_output_shape(
    stimulated_emission_factor, lines, level_number_density
):
    assert stimulated_emission_factor.shape == (
        len(lines),
        level_number_density.shape[1],
    )


def test_stimulated_emission_factor_lte_values_not_above_one(
    stimulated_emission_factor,
):
    assert np.all(stimulated_emission_factor <= 1.0)


def test_stimulated_emission_factor_regression(
    stimulated_emission_factor, regression_data
):
    expected = regression_data.sync_ndarray(
        stimulated_emission_factor
    )
    np.testing.assert_allclose(
        stimulated_emission_factor, expected, atol=0, rtol=1e-15
    )


def _two_level_inputs(n_lower, n_upper, g_lower, g_upper, metastable_upper=False):
    """Minimal 2-level, 1-zone, 1-line system for edge-case unit tests."""
    levels_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    lines_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0, 1)],
        names=["atomic_number", "ion_number", "level_number_lower", "level_number_upper"],
    )
    return dict(
        g=pd.Series([g_lower, g_upper], index=levels_index),
        level_number_density=pd.DataFrame(
            [[n_lower], [n_upper]], index=levels_index, columns=[0]
        ),
        metastability=pd.Series([False, metastable_upper], index=levels_index),
        lines=pd.DataFrame({"wavelength": [6562.8]}, index=lines_index),
        lines_lower_level_index=np.array([0], dtype=np.int64),
        lines_upper_level_index=np.array([1], dtype=np.int64),
    )


def test_stimulated_emission_factor_n_lower_zero_gives_zero():
    inputs = _two_level_inputs(n_lower=0.0, n_upper=1e10, g_lower=2.0, g_upper=4.0)
    actual = StimulatedEmissionFactor().calculate(**inputs)
    assert actual[0, 0] == 0.0


def test_stimulated_emission_factor_metastable_upper_clamps_negative_to_zero():
    # g_lower > g_upper forces a negative raw factor; metastable flag must clamp it
    inputs = _two_level_inputs(
        n_lower=1e10, n_upper=3e10, g_lower=4.0, g_upper=2.0, metastable_upper=True
    )
    actual = StimulatedEmissionFactor().calculate(**inputs)
    assert actual[0, 0] == 0.0


def test_stimulated_emission_factor_nlte_species_clamps_negative_to_zero():
    # Same inverted population; NLTE flag must clamp the negative factor
    inputs = _two_level_inputs(
        n_lower=1e10, n_upper=3e10, g_lower=4.0, g_upper=2.0
    )
    actual = StimulatedEmissionFactor(nlte_species={(1, 0)}).calculate(**inputs)
    assert actual[0, 0] == 0.0


