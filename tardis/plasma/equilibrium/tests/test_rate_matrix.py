import astropy.units as u
import numpy as np
import pandas.testing as pdt
import pytest

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rate_matrix import IonRateMatrix, RateMatrix
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


def test_rate_matrix_solver(
    new_chianti_atomic_dataset_si,
    rate_solver_list,
    collisional_simulation_state,
    regression_data,
):
    rate_matrix_solver = RateMatrix(
        rate_solver_list, new_chianti_atomic_dataset_si.levels
    )

    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0, collisional_simulation_state.t_radiative, 1e6 * u.g / u.cm**3
    )

    actual = rate_matrix_solver.solve(rad_field, electron_dist)

    expected = regression_data.sync_dataframe(actual)

    pdt.assert_frame_equal(actual, expected)


@pytest.mark.parametrize("charge_conservation", [True, False])
def test_ion_rate_matrix_solver(
    radiative_ionization_solver,
    collisional_ionization_solver,
    collisional_simulation_state,
    charge_conservation,
    regression_data,
):
    rate_matrix_solver = IonRateMatrix(
        radiative_ionization_solver, collisional_ionization_solver
    )

    rad_field = DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros_like(collisional_simulation_state.t_radiative),
    )
    electron_dist = ThermalElectronEnergyDistribution(
        0, collisional_simulation_state.t_radiative, 1e6 * u.g / u.cm**3
    )

    actual = rate_matrix_solver.solve(
        rad_field,
        electron_dist,
        lte_level_pop,
        level_pop,
        lte_ion_pop,
        ion_pop,
        charge_conservation,
    )

    expected = regression_data.sync_dataframe(actual)

    pdt.assert_frame_equal(actual, expected)
