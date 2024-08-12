import pytest
import pandas as pd
import numpy as np
from astropy import units as u

import numpy.testing as npt
import copy

from tardis.plasma.assembly.base import (
    PlasmaSolverFactory,
    convert_species_to_multi_index,
)
from tardis.plasma.radiation_field import planck_rad_field
from tardis.plasma.properties.atomic import YgData
from tardis.plasma.properties.plasma_input import ContinuumInteractionSpecies
from tardis.plasma.detailed_balance.rates import (
    RadiativeRatesSolver,
    UpsilonRegemorterSolver,
    UpsilonCMFGENSolver,
    ThermalCollisionalRateSolver,
)


@pytest.fixture
def all_legacy_cmfgen_yg_data(
    nlte_atomic_dataset,
):  # using christian's old implementation
    atom_data = copy.deepcopy(nlte_atomic_dataset)

    # almost all settings are irrelevant for collisional strength data
    number_densities = pd.DataFrame({2: [1]}).T
    temperatures = [10000] * u.K
    dilution_factor = np.array([1])
    time_explosion = 5 * u.day
    dilute_planck_rad_field = planck_rad_field.DilutePlanckianRadiationField(
        temperatures, dilution_factor
    )
    plasma_solver_factory = PlasmaSolverFactory(atom_data)
    species_mindex = convert_species_to_multi_index(["He I"])
    # plasma_solver_factory.continuum_interaction_species = ["He I"]
    plasma_solver_factory.line_interaction_type = "macroatom"
    plasma_solver_factory.prepare_factory([2])
    plasma_solver_factory.plasma_modules += [
        YgData,
        ContinuumInteractionSpecies,
    ]
    plasma_solver = plasma_solver_factory.assemble(
        number_densities,
        dilute_planck_rad_field,
        time_explosion,
        continuum_interaction_species=species_mindex,
    )
    available_yg_data = plasma_solver.yg_data.loc[
        atom_data.yg_data.loc[(2, 0, slice(None), slice(None)), :].index
    ]
    approximated_yg_data = plasma_solver.yg_data.loc[
        ~plasma_solver.yg_data.index.isin(atom_data.yg_data.index)
    ]
    return available_yg_data, approximated_yg_data


def test_legacy_cmfgen_collisional_strengths(
    all_legacy_cmfgen_yg_data, nlte_atomic_dataset, regression_data
):
    legacy_cmfgen_yg_data, approximated_cmfgen_yg_data = (
        all_legacy_cmfgen_yg_data
    )

    # This is testing againt the old setup
    radiative_transitions = nlte_atomic_dataset.lines.loc[
        (2, 0, slice(None), slice(None)), :
    ]

    collision_strengths_regemorter_solver = UpsilonRegemorterSolver(
        radiative_transitions.loc[approximated_cmfgen_yg_data.index]
    )
    collision_strengths_cmfgen_solver = UpsilonCMFGENSolver
    new_regemorter_collision_strengths = (
        collision_strengths_regemorter_solver.solve(
            t_electrons=legacy_cmfgen_yg_data.columns.values * u.K
        )
    )
    npt.assert_allclose(
        new_regemorter_collision_strengths.values, approximated_cmfgen_yg_data
    )  # residuals are ~1e-8 not sure if that is good enough
