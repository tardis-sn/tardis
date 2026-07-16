import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
from astropy import units as u

from tardis import constants as const
from tardis.iip_plasma.properties.continuum import (
    CollIonRateCoeff,
    CollRecombRateCoeff,
)
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rates.collisional_ionization_rates import (
    CollisionalIonizationRateSolver,
)
from tardis.plasma.equilibrium.rates.collisional_ionization_strengths import (
    CollisionalIonizationSeaton,
)


@pytest.fixture
def real_photoionization_data(nlte_atom_data):
    """Regression photoionization data spanning several ion charges."""
    return nlte_atom_data.photoionization_data.sort_values(
        ["atomic_number", "ion_number", "level_number", "nu"]
    )


@pytest.fixture
def electron_temperatures():
    return np.array([9000.0, 12000.0]) * u.K


@pytest.fixture
def level_to_ion_factor(
    real_photoionization_data, nlte_atom_data, electron_temperatures
):
    """LTE level-to-ion factors from the regression atomic data."""
    levels = nlte_atom_data.levels
    partition_functions = levels["g"].groupby(
        level=["atomic_number", "ion_number"]
    ).sum()
    electron_factor = (
        2
        * (
            2
            * np.pi
            * const.m_e.cgs.value
            * const.k_B.cgs.value
            * electron_temperatures.to_value(u.K)
            / const.h.cgs.value**2
        )
        ** 1.5
    )
    rows = []
    for atomic_number, ion_number, level_number in (
        real_photoionization_data.index.unique()
    ):
        level = levels.loc[(atomic_number, ion_number, level_number)]
        ionization_energy = nlte_atom_data.ionization_data.loc[
            (atomic_number, ion_number + 1)
        ]
        saha_factor = (
            2
            * electron_factor
            * partition_functions.loc[(atomic_number, ion_number + 1)]
            / partition_functions.loc[(atomic_number, ion_number)]
            * np.exp(
                -ionization_energy
                / (const.k_B.cgs.value * electron_temperatures.to_value(u.K))
            )
        )
        rows.append(
            level["g"]
            * np.exp(
                -level["energy"]
                / (const.k_B.cgs.value * electron_temperatures.to_value(u.K))
            )
            / (
                saha_factor
                * partition_functions.loc[(atomic_number, ion_number)]
            )
        )
    return pd.DataFrame(
        rows,
        index=real_photoionization_data.index.unique(),
    )


def test_seaton_thresholds_and_coefficients_match_analytic_expression(
    real_photoionization_data, electron_temperatures
):
    threshold_data = real_photoionization_data.groupby(
        level=["atomic_number", "ion_number", "level_number"]
    ).first()
    u0 = (
        threshold_data["nu"].to_numpy()[:, np.newaxis]
        * const.h.cgs.value
        / (const.k_B.cgs.value * electron_temperatures.to_value(u.K))
    )
    charge_factor = np.select(
        [
            threshold_data.index.get_level_values("ion_number") == 0,
            threshold_data.index.get_level_values("ion_number") == 1,
        ],
        [0.1, 0.2],
        default=0.3,
    )[:, np.newaxis]
    expected = (
        1.55e13
        * threshold_data["x_sect"].to_numpy()[:, np.newaxis]
        * charge_factor
        * np.exp(-u0)
        / u0
        / np.sqrt(electron_temperatures.to_value(u.K))
    )

    actual = CollisionalIonizationSeaton(real_photoionization_data).solve(
        electron_temperatures
    )
    npt.assert_allclose(actual.to_numpy(), expected, rtol=2e-14)
    assert np.all(actual.to_numpy() > 0)

    # The IIP property is a compatibility comparison; the expression above
    # is the independent Seaton/Hummer-Mihalas rate oracle.
    iip = CollIonRateCoeff(None).calculate(
        real_photoionization_data, electron_temperatures.to_value(u.K)
    )
    pd.testing.assert_frame_equal(actual, iip, check_names=False)


def test_seaton_temperature_dependence_matches_threshold_exponential(
    real_photoionization_data, electron_temperatures
):
    actual = CollisionalIonizationSeaton(real_photoionization_data).solve(
        electron_temperatures
    )
    threshold_data = real_photoionization_data.groupby(
        level=["atomic_number", "ion_number", "level_number"]
    ).first()
    u0 = (
        threshold_data["nu"].to_numpy()[:, np.newaxis]
        * const.h.cgs.value
        / const.k_B.cgs.value
        / electron_temperatures.to_value(u.K)
    )
    expected_ratio = (
        np.sqrt(electron_temperatures[0] / electron_temperatures[1])
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
    real_photoionization_data,
    electron_temperatures,
    level_to_ion_factor,
):
    partition_function = 1.0
    level_boltzmann_factor = pd.DataFrame(
        np.ones_like(level_to_ion_factor),
        index=level_to_ion_factor.index,
    )
    low_density = ThermalElectronEnergyDistribution(
        0 * u.erg,
        electron_temperatures,
        np.full(2, 1.0e9) / u.cm**3,
    )
    high_density = ThermalElectronEnergyDistribution(
        0 * u.erg,
        electron_temperatures,
        np.full(2, 2.0e9) / u.cm**3,
    )
    solver = CollisionalIonizationRateSolver(real_photoionization_data)
    ion_low, recomb_low = solver.solve(
        low_density,
        level_to_ion_factor,
        partition_function,
        level_boltzmann_factor,
    )
    ion_high, recomb_high = solver.solve(
        high_density,
        level_to_ion_factor,
        partition_function,
        level_boltzmann_factor,
    )

    npt.assert_allclose(ion_high.to_numpy(), 2 * ion_low.to_numpy())
    npt.assert_allclose(recomb_high.to_numpy(), 4 * recomb_low.to_numpy())
    assert np.all(ion_low.to_numpy() > 0)
    assert np.all(recomb_low.to_numpy() > 0)


def test_three_body_recombination_uses_lte_detailed_balance_factor(
    real_photoionization_data,
    electron_temperatures,
    level_to_ion_factor,
):
    coll_ion = CollisionalIonizationSeaton(real_photoionization_data).solve(
        electron_temperatures
    )
    expected = CollRecombRateCoeff(None).calculate(
        level_to_ion_factor, coll_ion
    )
    actual = coll_ion.multiply(level_to_ion_factor)

    npt.assert_allclose(actual.to_numpy(), expected.to_numpy(), rtol=1e-14)
    assert np.all(actual.to_numpy() > 0)
