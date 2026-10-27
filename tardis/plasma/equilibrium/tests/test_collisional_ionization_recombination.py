import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy import units as u
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis import constants as const
from tardis.io.atom_data import AtomData
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rates.collisional_ionization_rates import (
    CollisionalIonizationRateSolver,
)
from tardis.plasma.equilibrium.rates.collisional_ionization_strengths import (
    CollisionalIonizationSeaton,
)
from tardis.plasma.properties.general import BetaElectron, ThermalGElectron
from tardis.plasma.properties.ion_population import (
    SahaFactor,
    ThermalPhiSahaLTE,
)
from tardis.plasma.properties.partition_function import (
    ThermalLevelBoltzmannFactorLTE,
    ThermalLTEPartitionFunction,
)

REFERENCE_ELECTRON_TEMPERATURES = np.array([9000.0, 12000.0]) * u.K
REFERENCE_SEATON_RATE_COEFFICIENT = 1.55e13
REFERENCE_LOW_ELECTRON_DENSITY_CM3 = 1.0e9
REFERENCE_HIGH_ELECTRON_DENSITY_CM3 = 2.0e9


@pytest.fixture
def real_photoionization_data(nlte_atom_data: AtomData) -> pd.DataFrame:
    """Regression photoionization data spanning several ion charges."""
    return nlte_atom_data.photoionization_data.sort_values(
        ["atomic_number", "ion_number", "level_number", "nu"]
    )


@pytest.fixture
def level_to_ion_factor(
    real_photoionization_data: pd.DataFrame,
    nlte_atom_data: AtomData,
) -> pd.DataFrame:
    """LTE level-to-ion factors from the regression atomic data."""
    beta_electron = BetaElectron(None).calculate(
        REFERENCE_ELECTRON_TEMPERATURES.to_value(u.K)
    )
    level_boltzmann_factor = ThermalLevelBoltzmannFactorLTE(None).calculate(
        nlte_atom_data.levels["energy"],
        nlte_atom_data.levels["g"],
        beta_electron,
        nlte_atom_data.levels.index,
    )
    partition_function = ThermalLTEPartitionFunction(None).calculate(
        level_boltzmann_factor
    )
    thermal_phi_lte = ThermalPhiSahaLTE(None).calculate(
        ThermalGElectron(None).calculate(beta_electron),
        beta_electron,
        partition_function,
        nlte_atom_data.ionization_data,
    )
    return (
        SahaFactor(None)
        .calculate(thermal_phi_lte, level_boltzmann_factor, partition_function)
        .loc[real_photoionization_data.index.unique()]
    )


def test_seaton_thresholds_and_coefficients_match_analytic_expression(
    real_photoionization_data: pd.DataFrame,
    regression_data: RegressionData,
) -> None:
    actual = CollisionalIonizationSeaton(real_photoionization_data).solve(
        REFERENCE_ELECTRON_TEMPERATURES
    )
    expected = regression_data.sync_dataframe(actual, key="frame_0")
    pdt.assert_frame_equal(actual, expected, check_names=False)
    threshold_data = real_photoionization_data.groupby(
        level=["atomic_number", "ion_number", "level_number"]
    ).first()
    u0 = (
        threshold_data["nu"].to_numpy()[:, None]
        * const.h.cgs.value
        / (const.k_B.cgs.value * REFERENCE_ELECTRON_TEMPERATURES.to_value(u.K))
    )
    charge_factor = np.select(
        [
            threshold_data.index.get_level_values("ion_number") == 0,
            threshold_data.index.get_level_values("ion_number") == 1,
        ],
        [0.1, 0.2],
        default=0.3,
    )[:, None]
    # Hubeny & Mihalas Eq. 9.60 cgs prefactor (K**0.5 cm / s), p. 276:
    # https://books.google.com/books?id=VA_rAwAAQBAJ&pg=PA276
    analytic = (
        REFERENCE_SEATON_RATE_COEFFICIENT
        * threshold_data["x_sect"].to_numpy()[:, None]
        * charge_factor
        * np.exp(-u0)
        / u0
        / np.sqrt(REFERENCE_ELECTRON_TEMPERATURES.to_value(u.K))
    )
    npt.assert_allclose(actual.to_numpy(), analytic, rtol=2e-14)


def test_seaton_temperature_dependence_matches_threshold_exponential(
    real_photoionization_data: pd.DataFrame,
) -> None:
    actual = CollisionalIonizationSeaton(real_photoionization_data).solve(
        REFERENCE_ELECTRON_TEMPERATURES
    )
    threshold_data = real_photoionization_data.groupby(
        level=["atomic_number", "ion_number", "level_number"]
    ).first()
    u0 = (
        threshold_data["nu"].to_numpy()[:, None]
        * const.h.cgs.value
        / const.k_B.cgs.value
        / REFERENCE_ELECTRON_TEMPERATURES.to_value(u.K)
    )
    expected_ratio = (
        np.sqrt(
            REFERENCE_ELECTRON_TEMPERATURES[0]
            / REFERENCE_ELECTRON_TEMPERATURES[1]
        )
        * np.exp(u0[:, 0] - u0[:, 1])
        * u0[:, 0]
        / u0[:, 1]
    )
    npt.assert_allclose(
        (actual.iloc[:, 1] / actual.iloc[:, 0]).to_numpy(),
        expected_ratio,
        rtol=2e-14,
    )


def test_collisional_rates_scale_with_electron_density(
    real_photoionization_data: pd.DataFrame,
    level_to_ion_factor: pd.DataFrame,
) -> None:
    partition_function = pd.DataFrame(
        1.0,
        index=level_to_ion_factor.index,
        columns=level_to_ion_factor.columns,
    )
    level_boltzmann_factor = pd.DataFrame(
        np.ones_like(level_to_ion_factor),
        index=level_to_ion_factor.index,
    )
    solver = CollisionalIonizationRateSolver(real_photoionization_data)
    rates = []
    for density in (
        REFERENCE_LOW_ELECTRON_DENSITY_CM3,
        REFERENCE_HIGH_ELECTRON_DENSITY_CM3,
    ):
        distribution = ThermalElectronEnergyDistribution(
            0 * u.erg,
            REFERENCE_ELECTRON_TEMPERATURES,
            np.full(2, density) / u.cm**3,
        )
        rates.append(
            solver.solve(
                distribution,
                level_to_ion_factor,
                partition_function,
                level_boltzmann_factor,
            )
        )
    npt.assert_allclose(rates[1][0].to_numpy(), 2 * rates[0][0].to_numpy())
    npt.assert_allclose(rates[1][1].to_numpy(), 4 * rates[0][1].to_numpy())


def test_three_body_recombination_uses_lte_detailed_balance_factor(
    real_photoionization_data: pd.DataFrame,
    level_to_ion_factor: pd.DataFrame,
    regression_data: RegressionData,
) -> None:
    actual = (
        CollisionalIonizationSeaton(real_photoionization_data)
        .solve(REFERENCE_ELECTRON_TEMPERATURES)
        .multiply(level_to_ion_factor)
    )
    expected = regression_data.sync_dataframe(
        pd.DataFrame(actual.to_numpy()), key="allclose_0"
    )
    # The snapshot is the legacy IIP detailed-balance calculation.
    npt.assert_allclose(actual.to_numpy(), expected.to_numpy(), rtol=1e-14)
    assert np.all(actual.to_numpy() > 0)
