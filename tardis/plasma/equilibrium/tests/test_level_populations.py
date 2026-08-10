import astropy.units as u
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.io.atom_data import AtomData
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.ion_populations import (
    AnalyticIonPopulationSolver,
)
from tardis.plasma.equilibrium.level_populations import LevelPopulationSolver
from tardis.plasma.equilibrium.rate_matrix import (
    IonRateMatrix,
    LevelRateMatrix,
)
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
    RadiativeRatesSolver,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.properties.atomic import IonizationData, Levels
from tardis.plasma.properties.general import BetaRadiation, GElectron
from tardis.plasma.properties.ion_population import IonNumberDensity, PhiSahaLTE
from tardis.plasma.properties.level_population import LevelNumberDensity
from tardis.plasma.properties.partition_function import (
    LevelBoltzmannFactorLTE,
    PartitionFunction,
)
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


class TestLevelPopulationSolver:
    @pytest.fixture(autouse=True)
    def setup(
        self,
        radiative_rate_solver,
        collisional_rate_solver,
        new_chianti_atomic_dataset_si,
        collisional_simulation_state,
    ):
        rad_field = DilutePlanckianRadiationField(
            collisional_simulation_state.t_radiative,
            dilution_factor=np.zeros_like(
                collisional_simulation_state.t_radiative
            ),
        )
        electron_dist = ThermalElectronEnergyDistribution(
            0, collisional_simulation_state.t_radiative, 1e6 * u.g / u.cm**3
        )

        radiative_rates = radiative_rate_solver.solve(rad_field)
        collisional_rates = collisional_rate_solver.solve(
            electron_dist.temperature
        )
        rates_matrices = LevelRateMatrix(
            new_chianti_atomic_dataset_si.levels
        ).solve(
            radiative_rates,
            collisional_rates,
            electron_dist.number_density.value,
        )
        self.solver = LevelPopulationSolver(
            rates_matrices, new_chianti_atomic_dataset_si.levels
        )

    def test_solve(self, regression_data):
        """Test the solve method."""
        result = self.solver.solve()
        expected_populations = regression_data.sync_dataframe(result)
        pdt.assert_frame_equal(result, expected_populations, atol=0, rtol=1e-15)

        for species_id in self.solver.rates_matrices.index:
            species_matrices = self.solver.rates_matrices.loc[species_id]
            species_populations = result.loc[species_id]
            for shell in result.columns:
                matrix = species_matrices[shell]
                population = species_populations[shell].to_numpy()
                balance = np.zeros(matrix.shape[0])
                balance[0] = 1.0

                np.testing.assert_allclose(
                    matrix @ population,
                    balance,
                    rtol=1e-12,
                    atol=1e-14,
                )
                np.testing.assert_allclose(
                    population.sum(),
                    1.0,
                    rtol=1e-12,
                    atol=1e-14,
                )
                assert np.all(population >= 0.0)
                assert np.isfinite(np.linalg.cond(matrix))


def test_equilibrium_rate_matrices_converge_to_equilibrium_lte(
    nlte_atomic_dataset: AtomData,
) -> None:
    """Verify all equilibrium rate matrices recover their own LTE limit."""
    atomic_data = nlte_atomic_dataset
    species = (1, 0)
    selected_atoms = pd.Index([species[0]], name="atomic_number")
    temperatures = np.array([10000.0, 20000.0]) * u.K
    # Use the collision-dominated limit while keeping ionization and
    # recombination active in the ion rate matrix.
    elemental_number_density = pd.DataFrame(
        1.0e25,
        index=selected_atoms,
        columns=range(len(temperatures)),
    )

    radiation_field = DilutePlanckianRadiationField(
        temperatures,
        dilution_factor=np.ones(len(temperatures)),
    )
    levels, excitation_energy, _, statistical_weights = Levels(None).calculate(
        atomic_data,
        selected_atoms,
    )
    beta_radiation = BetaRadiation(None).calculate(temperatures.to_value(u.K))
    lte_boltzmann_factors = LevelBoltzmannFactorLTE(None).calculate(
        excitation_energy,
        statistical_weights,
        beta_radiation,
        levels,
    )
    lte_partition_function = PartitionFunction(None).calculate(
        lte_boltzmann_factors
    )
    ionization_data = IonizationData(None).calculate(
        atomic_data,
        selected_atoms,
    )
    saha_factor = PhiSahaLTE(None).calculate(
        GElectron(None).calculate(beta_radiation),
        beta_radiation,
        lte_partition_function,
        ionization_data,
    )
    hydrogen_saha_factor = saha_factor.loc[
        (species[0], species[1] + 1)
    ].to_numpy()
    hydrogen_density = elemental_number_density.loc[species[0]].to_numpy()
    electron_densities = pd.Series(
        2.0
        * hydrogen_density
        / (1.0 + np.sqrt(1.0 + 4.0 * hydrogen_density / hydrogen_saha_factor))
    )
    lte_ion_populations, _ = IonNumberDensity(
        None,
        electron_densities=electron_densities,
    ).calculate(saha_factor, lte_partition_function, elemental_number_density)
    lte_electron_densities = (
        lte_ion_populations
        * lte_ion_populations.index.get_level_values("ion_number").to_numpy()[
            :, None
        ]
    ).sum()
    lte_level_populations = LevelNumberDensity(None).calculate(
        lte_boltzmann_factors,
        lte_ion_populations,
        levels,
        lte_partition_function,
    )

    radiative_transitions = atomic_data.lines.loc[
        (species[0], species[1], slice(None), slice(None)),
        :,
    ]
    collision_strengths = atomic_data.yg_data.loc[
        (species[0], species[1], slice(None), slice(None)),
        :,
    ]
    radiative_rate_solver = RadiativeRatesSolver(radiative_transitions)
    collisional_rate_solver = ThermalCollisionalRateSolver(
        atomic_data.levels,
        radiative_transitions,
        atomic_data.collision_data_temperatures,
        collision_strengths,
        collision_strengths_type="cmfgen",
    )
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg,
        temperatures,
        electron_densities.to_numpy() / u.cm**3,
    )
    level_matrices = LevelRateMatrix(atomic_data.levels).solve(
        radiative_rate_solver.solve(radiation_field),
        collisional_rate_solver.solve(electron_distribution.temperature),
        electron_distribution.number_density.value,
    )
    nlte_level_population_fractions = (
        LevelPopulationSolver(
            level_matrices,
            atomic_data.levels,
        )
        .solve()
        .loc[species]
    )
    lte_level_population_fractions = lte_level_populations.loc[species].divide(
        lte_ion_populations.loc[species],
        axis=1,
    )

    np.testing.assert_allclose(
        nlte_level_population_fractions.to_numpy(),
        lte_level_population_fractions.to_numpy(),
        # Independent rate and Boltzmann calculations agree to about 5e-8.
        # Increased to deal with Mac vs Linux numerical differences.
        rtol=5e-7,
        atol=0.0,
    )

    estimated_level_populations = lte_level_populations.copy()
    estimated_level_populations.loc[species] = (
        nlte_level_population_fractions.multiply(
            lte_ion_populations.loc[species],
            axis=1,
        ).to_numpy()
    )
    photoionization_data = atomic_data.photoionization_data.loc[
        (species[0], species[1], slice(None)),
        :,
    ].sort_values("nu")
    ion_population_solver = AnalyticIonPopulationSolver(
        IonRateMatrix(),
        AnalyticPhotoionizationRateSolver(photoionization_data),
        CollisionalIonizationRateSolver(photoionization_data),
        elemental_number_density,
        radiation_field=radiation_field,
        thermal_electron_energy_distribution=electron_distribution,
        lte_level_population=lte_level_populations,
        estimated_level_population=estimated_level_populations,
        lte_ion_population=lte_ion_populations,
        estimated_ion_population=lte_ion_populations.copy(),
        partition_function=lte_partition_function,
        boltzmann_factor=lte_boltzmann_factors,
    )
    nlte_ion_populations, nlte_electron_densities = (
        ion_population_solver.solve()
    )

    np.testing.assert_allclose(
        nlte_ion_populations.to_numpy(),
        lte_ion_populations.to_numpy(),
        rtol=1e-8,
        atol=0.0,
    )
    np.testing.assert_allclose(
        nlte_electron_densities.to_numpy(),
        lte_electron_densities.to_numpy(),
        rtol=1e-8,
        atol=0.0,
    )


def test_level_population_solver_solves_raw_level_block() -> None:
    levels = pd.DataFrame(
        {"energy": [0.0, 1.0]},
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
    )
    raw_rate_matrices = pd.DataFrame(
        [[np.array([[-2.0, 3.0], [2.0, -3.0]])]],
        index=pd.MultiIndex.from_tuples(
            [(1, 0)], names=["atomic_number", "ion_number"]
        ),
        columns=[0],
    )

    solver = LevelPopulationSolver(raw_rate_matrices, levels)
    populations = solver.solve_raw(raw_rate_matrices)

    np.testing.assert_allclose(populations.loc[(1, 0), 0], [0.6, 0.4])
