import copy

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy import units as u

from tardis.io.atom_data import AtomData
from tardis.plasma.assembly.base import (
    PlasmaSolverFactory,
    convert_species_to_multi_index,
)
from tardis.plasma.equilibrium.rates import (
    #    UpsilonCMFGENSolver,
    ThermalCollisionalRateSolver,
    #    RadiativeRatesSolver,
    UpsilonRegemorterSolver,
)
from tardis.plasma.properties.atomic import YgData, YgInterpolator
from tardis.plasma.properties.continuum_processes import (
    CollDeexcRateCoeff,
    CollExcRateCoeff,
)
from tardis.plasma.properties.general import BetaElectron
from tardis.plasma.properties.partition_function import (
    ThermalLevelBoltzmannFactorLTE,
)
from tardis.plasma.properties.plasma_input import ContinuumInteractionSpecies
from tardis.plasma.radiation_field import planck_rad_field


@pytest.fixture
def legacy_cmfgen_collision_rate_plasma_solver(nlte_atomic_dataset):
    atom_data = copy.deepcopy(nlte_atomic_dataset)
    # almost all settings are irrelevant for collisional strength data
    number_densities = pd.DataFrame({1: [1, 1]}).T
    temperatures = [10000, 20000] * u.K
    dilution_factor = np.array([1, 1])
    time_explosion = 5 * u.day
    dilute_planck_rad_field = planck_rad_field.DilutePlanckianRadiationField(
        temperatures, dilution_factor
    )
    plasma_solver_factory = PlasmaSolverFactory(atom_data)

    # plasma_solver_factory.continuum_interaction_species = ["He I"]
    plasma_solver_factory.line_interaction_type = "macroatom"
    plasma_solver_factory.prepare_factory(
        [1], "tardis.plasma.properties.legacy_property_collections"
    )
    plasma_solver_factory.plasma_modules += [
        YgData,
        ContinuumInteractionSpecies,
        CollExcRateCoeff,
        CollDeexcRateCoeff,
        YgInterpolator,
        ThermalLevelBoltzmannFactorLTE,
        BetaElectron,
    ]
    species_mindex = convert_species_to_multi_index(["H I"])
    return plasma_solver_factory.assemble(
        number_densities,
        dilute_planck_rad_field,
        time_explosion,
        continuum_interaction_species=species_mindex,
    )


@pytest.fixture
def new_chianti_atomic_dataset(tardis_regression_path):
    atomic_data_fname = (
        tardis_regression_path / "atom_data" / "new_kurucz_cd23_chianti_H_He.h5"
    )
    return AtomData.from_hdf(atomic_data_fname)


@pytest.fixture
def legacy_chianti_collision_rate_plasma_solver(atomic_dataset):
    atom_data = copy.deepcopy(atomic_dataset)
    atom_data.prepare_atom_data([1], "macroatom", [(1, 0)], [])
    return atom_data.nlte_data.get_collision_matrix(
        (1, 0), np.array([10000, 20000])
    )


def test_legacy_cmfgen_collisional_strengths(
    legacy_cmfgen_collision_rate_plasma_solver,
    nlte_atomic_dataset,
    regression_data,
):
    # using christian's old implementation
    plasma_solver = legacy_cmfgen_collision_rate_plasma_solver
    atom_data = copy.deepcopy(nlte_atomic_dataset)
    legacy_cmfgen_yg_data = plasma_solver.yg_data.loc[
        atom_data.yg_data.loc[(1, 0, slice(None), slice(None)), :].index
    ]
    approximated_cmfgen_yg_data = plasma_solver.yg_data.loc[
        ~plasma_solver.yg_data.index.isin(atom_data.yg_data.index)
    ]

    # This is testing againt the old setup
    radiative_transitions = atom_data.lines.loc[
        (1, 0, slice(None), slice(None)), :
    ]

    collision_strengths_regemorter_solver = UpsilonRegemorterSolver(
        radiative_transitions.loc[approximated_cmfgen_yg_data.index]
    )

    new_regemorter_collision_strengths = (
        collision_strengths_regemorter_solver.solve(
            t_electrons=legacy_cmfgen_yg_data.columns.values * u.K
        )
    )
    npt.assert_allclose(
        new_regemorter_collision_strengths.values,
        approximated_cmfgen_yg_data,
    )  # residuals are ~1e-8 not sure if that is good enough
    # Not comparing to the yg_data as they are saved differently


def test_thermal_collision_rates(
    legacy_cmfgen_collision_rate_plasma_solver,
    nlte_atomic_dataset,
    regression_data,
):
    atom_data = copy.deepcopy(nlte_atomic_dataset)
    radiative_transitions = atom_data.lines.loc[
        (1, 0, slice(None), slice(None)), :
    ]

    collision_strengths = atom_data.yg_data.loc[
        (1, 0, slice(None), slice(None)), :
    ]
    collision_strengths_temperatures = atom_data.collision_data_temperatures

    therm_coll_rate_solver = ThermalCollisionalRateSolver(
        atom_data.levels,
        radiative_transitions,
        collision_strengths_temperatures,
        collision_strengths,
        collision_strengths_type="cmfgen",
        collisional_strength_approximation="regemorter",
    )
    coll_rates_coeff = therm_coll_rate_solver.solve([10000, 20000] * u.K)
    pdt.assert_frame_equal(
        coll_rates_coeff.iloc[:435],
        legacy_cmfgen_collision_rate_plasma_solver.coll_exc_coeff,
        check_names=False,
    )
    pdt.assert_frame_equal(
        coll_rates_coeff.iloc[435:],
        legacy_cmfgen_collision_rate_plasma_solver.coll_deexc_coeff.swaplevel(
            "level_number_lower", "level_number_upper"
        ),
        check_names=False,
    )


# Add chianti tests
def test_legacy_chianti_collisional_strengths(
    legacy_chianti_collision_rate_plasma_solver,
    new_chianti_atomic_dataset,
    regression_data,
):
    collision_strengths = legacy_chianti_collision_rate_plasma_solver
    atom_data = copy.deepcopy(new_chianti_atomic_dataset)

    temperature = np.array([10000, 20000]) * u.K

    col_strengths = atom_data.collision_data.loc[
        (1, 0, slice(None), slice(None)), :
    ]
    radiative_transitions = atom_data.lines.loc[
        (1, 0, slice(None), slice(None)), :
    ]
    collisional_rate_solver = ThermalCollisionalRateSolver(
        atom_data.levels,
        radiative_transitions,
        temperature,
        col_strengths,
        "chianti",
    )
    chianti_collisional_rates = collisional_rate_solver.solve(temperature)

    npt.assert_allclose(
        collision_strengths[0, 1, :],
        chianti_collisional_rates.loc[1, 0, 0, 0, 1, 0],
        rtol=1e-4,
        atol=1e-13,
    )
