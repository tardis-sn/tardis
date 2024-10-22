import astropy.units as u
import numpy as np
import pytest

from tardis.plasma.radiation_field.planck_rad_field import (
    DilutePlanckianRadiationField,
)
from tardis.simulation.radiation_field_convergence import (
    RadiationFieldConvergenceSolver,
)


@pytest.fixture(scope="function")
def radiation_field():
    temperature = [1000, 2000, 3000] * u.K
    dilution_factor = np.array([0.1, 0.2, 0.3])
    radiation_field = DilutePlanckianRadiationField(
        temperature, dilution_factor
    )
    return radiation_field


@pytest.fixture(scope="function")
def partially_converged_radiation_field():
    temperature = [1000, 2000, 3000] * u.K
    dilution_factor = np.array([0.15, 0.25, 0.35])
    radiation_field = DilutePlanckianRadiationField(
        temperature, dilution_factor
    )
    return radiation_field


@pytest.mark.parametrize(
    "converge_separately, estimated_radiation_field_fixture, expected",
    [
        (True, "radiation_field", (True, True)),
        (True, "partially_converged_radiation_field", (True, False)),
        (False, "partially_converged_radiation_field", False),
    ],
)
def test_get_convergence_status(
    converge_separately,
    estimated_radiation_field_fixture,
    expected,
    request,
    strategy,
    radiation_field,
):
    estimated_radiation_field = request.getfixturevalue(
        estimated_radiation_field_fixture
    )
    solver = RadiationFieldConvergenceSolver(strategy)
    solver.converge_separately = converge_separately

    converged_result = solver.get_convergence_status(
        radiation_field, estimated_radiation_field, 3
    )
    assert converged_result == expected
