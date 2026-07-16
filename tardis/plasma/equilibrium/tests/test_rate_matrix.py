import astropy.units as u
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rate_matrix import (
    EquilibriumIonRateMatrix,
    LevelRateMatrix,
)
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


def test_rate_matrix_solver(
    new_chianti_atomic_dataset_si,
    rate_solver_list,
    collisional_simulation_state,
    regression_data,
):
    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0, collisional_simulation_state.t_radiative, 1e6 * u.g / u.cm**3
    )

    radiative_rates = rate_solver_list[0][0].solve(rad_field)
    collisional_rates = rate_solver_list[1][0].solve(electron_dist.temperature)
    actual = LevelRateMatrix(new_chianti_atomic_dataset_si.levels).solve(
        radiative_rates,
        collisional_rates,
        electron_dist.number_density.value,
    )

    expected = regression_data.sync_dataframe(actual)

    pdt.assert_frame_equal(actual, expected, atol=0, rtol=1e-15)


@pytest.mark.parametrize("charge_conservation", [True, False])
def test_ion_rate_matrix_solver(
    photoionization_rate_solver,
    collisional_ionization_rate_solver,
    collisional_simulation_state,
    mock_boltzmann_factor,
    charge_conservation,
    regression_data,
):
    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0, collisional_simulation_state.t_radiative, 1e6 * u.g / u.cm**3
    )

    lte_level_population = pd.DataFrame(
        data=np.ones((2, 20)) * 1e5,
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
    )

    lte_ion_population = pd.DataFrame(
        data=np.ones((2, 20)) * 1e5,
        index=pd.MultiIndex.from_tuples(
            [(1, 0), (1, 1)],
            names=["atomic_number", "ion_number"],
        ),
    )

    level_population = lte_level_population.copy() * 1.4
    ion_population = lte_ion_population.copy() * 3.0

    photoion_rates, recombination_rates = photoionization_rate_solver.solve(
        rad_field,
        electron_dist,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
        1.0,
        mock_boltzmann_factor,
    )
    level_to_ion_population_factor = lte_level_population / (
        lte_ion_population.values * electron_dist.number_density.value
    )
    collisional_ionization_rates, collisional_recombination_rates = (
        collisional_ionization_rate_solver.solve(
            electron_dist,
            level_to_ion_population_factor,
            1.0,
            mock_boltzmann_factor,
        )
    )
    actual = EquilibriumIonRateMatrix().solve(
        photoion_rates,
        recombination_rates,
        collisional_ionization_rates,
        collisional_recombination_rates,
        charge_conservation,
    )

    expected = regression_data.sync_dataframe(actual)

    pdt.assert_frame_equal(actual, expected, atol=0, rtol=1e-15)
