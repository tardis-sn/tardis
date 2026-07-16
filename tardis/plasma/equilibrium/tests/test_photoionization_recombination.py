import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
from astropy import units as u

from tardis import constants as const
from tardis.iip_plasma.properties.continuum import (
    PhotoIonRateCoeff,
    SpontRecombRateCoeff,
    StimRecombRateCoeff,
)
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
)
from tardis.plasma.equilibrium.rates.photoionization_strengths import (
    AnalyticPhotoionizationCoeffSolver,
    EstimatedPhotoionizationCoeffSolver,
    SpontaneousRecombinationCoeffSolver,
)
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.transport.montecarlo.estimators import init_estimators_continuum


@pytest.fixture
def photoionization_data(nlte_atom_data):
    return nlte_atom_data.photoionization_data.loc[(1, 0, [1])].sort_values(
        ["atomic_number", "ion_number", "level_number", "nu"]
    )


@pytest.fixture
def lyman_photoionization_data(nlte_atom_data):
    return nlte_atom_data.photoionization_data.loc[(1, 0, [0, 1])].sort_values(
        ["atomic_number", "ion_number", "level_number", "nu"]
    )


@pytest.fixture
def radiation_field():
    return DilutePlanckianRadiationField(
        np.array([10000.0, 12000.0]) * u.K, np.array([0.4, 0.8])
    )


def _level_frame(index, values):
    return pd.DataFrame(values, index=index, columns=[0, 1])


def test_analytic_photoionization_rates_match_iip(
    photoionization_data, radiation_field,
):
    photo_data = photoionization_data
    electron_temperature = np.array([9000.0, 11000.0]) * u.K
    standard_solver = AnalyticPhotoionizationCoeffSolver(photo_data)
    gamma, alpha_stim = standard_solver.solve(
        radiation_field, electron_temperature
    )

    iip_j_nu = PhotoIonRateCoeff._calculate_j_nus(
        photo_data, radiation_field.dilution_factor, radiation_field.temperature_kelvin
    )
    standard_j_nu = standard_solver.calculate_mean_intensity_photoionization_df(
        radiation_field
    )
    npt.assert_allclose(standard_j_nu.to_numpy(), iip_j_nu.to_numpy(), rtol=1e-12)

    iip_gamma = PhotoIonRateCoeff(None).calculate_from_radiation_field_model(
        photo_data,
        radiation_field.dilution_factor,
        radiation_field.temperature_kelvin,
    )
    iip_alpha_stim = StimRecombRateCoeff(
        None
    ).calculate_from_radiation_field_model(
        photo_data,
        radiation_field.dilution_factor,
        radiation_field.temperature_kelvin,
        None,
        electron_temperature.value,
        _level_frame(photo_data.index.unique(), np.ones((1, 2))),
    )
    # Provenance: gamma and alpha_stim are Lucy (2003) Eqs. 16 and 15,
    # respectively. These calls exercise the standard implementation in
    # photoionization_strengths.py against the corresponding IIP properties;
    # the test is intended to catch differences in grouping or normalization.
    # The implementations use the same integrand but different numerical
    # integration backends and constant sources.
    npt.assert_allclose(gamma.to_numpy(), iip_gamma.to_numpy(), rtol=2e-7)
    # The stimulated-recombination paths use different cgs constant sources
    # and quadrature implementations.
    npt.assert_allclose(alpha_stim.to_numpy(), iip_alpha_stim.to_numpy(), rtol=2e-6)


def test_zero_radiation_gives_zero_photoionization_and_stimulated_recombination(
    photoionization_data,
):
    photo_data = photoionization_data
    zero_field = DilutePlanckianRadiationField(
        np.array([10000.0, 12000.0]) * u.K, np.zeros(2)
    )
    gamma, alpha_stim = AnalyticPhotoionizationCoeffSolver(photo_data).solve(
        zero_field, np.array([9000.0, 11000.0]) * u.K
    )
    assert np.all(gamma.to_numpy() == 0.0)
    assert np.all(alpha_stim.to_numpy() == 0.0)


def test_spontaneous_recombination_is_positive_and_lyman_suppression_is_explicit(
    lyman_photoionization_data,
):
    photo_data = lyman_photoionization_data
    temperatures = np.array([10000.0, 12000.0]) * u.K
    standard_alpha = SpontaneousRecombinationCoeffSolver(photo_data).solve(
        temperatures
    )
    assert (standard_alpha.loc[(1, 0, 1)] >= 0).all()
    assert np.all(standard_alpha.loc[(1, 0, 0)] == 0.0)

    phi_lucy = _level_frame(photo_data.index.unique(), np.ones((2, 2)))
    iip_alpha = SpontRecombRateCoeff(
        type("IterationState", (), {"niter": 2, "niter_ly": 1})()
    ).calculate(photo_data, temperatures.value, phi_lucy)
    assert np.all(iip_alpha.loc[(1, 0, 0)] == 0.0)
    npt.assert_allclose(
        standard_alpha.loc[(1, 0, 1)].to_numpy(),
        iip_alpha.loc[(1, 0, 1)].to_numpy(),
        rtol=2e-4,
    )


def test_estimator_coefficients_reproduce_regression_inputs(tardis_regression_path):
    edge_index = pd.MultiIndex.from_tuples(
        [(1, 0, 1), (1, 0, 2)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    photo_ion_estimator = pd.read_hdf(
        tardis_regression_path
        / "tardis"
        / "workflows"
        / "tests"
        / "ctardis_photo_ion_estimator_after_mc.h5",
        key="data",
    ).loc[edge_index, [0, 1]]
    stim_recomb_estimator = pd.read_hdf(
        tardis_regression_path
        / "tardis"
        / "workflows"
        / "tests"
        / "ctardis_stim_recomb_estimator_after_mc.h5",
        key="data",
    ).loc[edge_index, [0, 1]]
    estimators = init_estimators_continuum(
        n_levels_bf_species_by_n_cells_tuple=(2, 2), n_cells=2
    )
    estimators.photo_ion_estimator[:] = photo_ion_estimator.to_numpy()
    estimators.stim_recomb_estimator[:] = stim_recomb_estimator.to_numpy()
    time_simulation = 2.0e5 * u.s
    volume = 3.0e30 * u.cm**3
    normalization = 1.0 / (
        time_simulation.to_value(u.s)
        * volume.to_value(u.cm**3)
        * const.h.cgs.value
    )
    gamma, alpha_stim = EstimatedPhotoionizationCoeffSolver(
        pd.Series([0, 1], index=edge_index)
    ).solve(estimators, time_simulation, volume)
    npt.assert_allclose(
        gamma.to_numpy(),
        photo_ion_estimator.to_numpy() * normalization,
    )
    npt.assert_allclose(
        alpha_stim.to_numpy(),
        stim_recomb_estimator.to_numpy() * normalization,
    )
    assert gamma.index.equals(edge_index)


def test_bound_free_heating_and_cooling_match_independent_quadrature(
    lyman_photoionization_data, radiation_field,
):
    photo_data = lyman_photoionization_data
    temperatures = np.array([10000.0, 12000.0])
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg, temperatures * u.K, np.ones(2) * 1.0e9 / u.cm**3
    )
    level_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    level_population = _level_frame(level_index, np.ones((2, 2)))
    ion_population = _level_frame(
        pd.MultiIndex.from_tuples([(1, 1)], names=["atomic_number", "ion_number"]),
        np.ones((1, 2)),
    )
    level_population_ratio = _level_frame(level_index, np.ones((2, 2)))
    heating, cooling = BoundFreeThermalRates(photo_data).solve(
        level_population,
        ion_population,
        electron_distribution,
        level_population_ratio,
        radiation_field,
    )
    thresholds = {
        level_number: photo_data.loc[(1, 0, level_number), "nu"].min()
        for level_number in [0, 1]
    }
    # Provenance: the heating and cooling integrands are Lucy (2003) Eqs. 58
    # and 59, as implemented by BoundFreeThermalRates and the matching IIP
    # continuum properties. Reconstruct them here directly from the real
    # regression cross-section grid so the test remains independent of the
    # solver's block-integration implementation.
    expected_heating = []
    expected_cooling = []
    for t_rad, t_electron, dilution in zip(
        [10000.0, 12000.0], temperatures, [0.4, 0.8]
    ):
        heating_integral = 0.0
        for level_number in [0, 1]:
            level_data = photo_data.loc[(1, 0, level_number)].sort_values("nu")
            frequencies = level_data["nu"].to_numpy()
            cross_section = level_data["x_sect"].to_numpy()
            mean_intensity = dilution * (
                2
                * const.h.cgs.value
                * frequencies**3
                / const.c.cgs.value**2
                * np.exp(
                    -const.h.cgs.value
                    * frequencies
                    / (const.k_B.cgs.value * t_rad)
                )
                / -np.expm1(
                    -const.h.cgs.value
                    * frequencies
                    / (const.k_B.cgs.value * t_rad)
                )
            )
            heating_integrand = (
                4
                * np.pi
                * cross_section
                * frequencies**3
                * const.h.cgs.value
                / const.c.cgs.value**2
                * (1 - thresholds[level_number] / frequencies)
                * mean_intensity
            )
            heating_integral += np.trapezoid(heating_integrand, frequencies)
        cooling_integrand = (
            8
            * np.pi
            * photo_data.loc[(1, 0, 1), "x_sect"].to_numpy()
            * frequencies**3
            * const.h.cgs.value
            / const.c.cgs.value**2
            * (1 - thresholds[1] / frequencies)
            * np.exp(-const.h.cgs.value * frequencies / (const.k_B.cgs.value * t_electron))
            * 1.0e9
        )
        expected_heating.append(heating_integral)
        expected_cooling.append(np.trapezoid(cooling_integrand, frequencies))
    npt.assert_allclose(heating.to_numpy(), expected_heating, rtol=1e-12)
    npt.assert_allclose(cooling.to_numpy(), expected_cooling, rtol=1e-12)
