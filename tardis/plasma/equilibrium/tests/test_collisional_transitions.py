import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy import units as u

from tardis import constants as const
from tardis.iip_plasma.properties.continuum import (
    CollDeexcRateCoeff,
    CollExcRateCoeff,
    Yg,
    YgInterpolator,
)
from tardis.iip_plasma.properties.partition_function import (
    LevelBoltzmannFactorLTETe,
)
from tardis.plasma.equilibrium.rates import (
    ThermalCollisionalRateSolver,
)


def _with_ion_source_and_destination_index(
    rate_coefficients: pd.DataFrame,
) -> pd.DataFrame:
    """Match the equilibrium solver's collisional-rate index layout."""
    rate_coefficients = rate_coefficients.copy()
    rate_coefficients.index.names = [
        "atomic_number",
        "ion_number",
        "level_number_source",
        "level_number_destination",
    ]
    rate_coefficients = rate_coefficients.reset_index()
    rate_coefficients["ion_number_source"] = rate_coefficients["ion_number"]
    rate_coefficients["ion_number_destination"] = rate_coefficients[
        "ion_number"
    ]
    return rate_coefficients.set_index(
        [
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
            "level_number_source",
            "level_number_destination",
        ]
    )


@pytest.fixture
def iip_collision_rate_coefficients(iip_atom_data: object) -> pd.DataFrame:
    """Calculate excitation and de-excitation rates through IIP properties."""
    temperatures_electron = np.array([10000.0, 20000.0])
    yg_interpolator = YgInterpolator(None)
    yg_interp, yg_allowed_index, yg_forbidden_index = (
        yg_interpolator.calculate(
            iip_atom_data.yg_data,
            iip_atom_data.collision_data_temperatures,
            iip_atom_data.lines.index,
        )
    )
    yg = Yg(None).calculate(
        yg_interp, iip_atom_data.yg_data.index, temperatures_electron
    )
    coll_exc_coeff, _ = CollExcRateCoeff(None).calculate(
        iip_atom_data.lines,
        iip_atom_data.levels.energy,
        temperatures_electron,
        yg,
        yg_allowed_index,
        yg_forbidden_index,
    )
    beta_electron = 1.0 / (const.k_B.cgs.value * temperatures_electron)
    lte_level_boltzmann_factor = LevelBoltzmannFactorLTETe.calculate(
        iip_atom_data.levels.energy,
        iip_atom_data.levels.g,
        beta_electron,
        iip_atom_data.levels.index,
    )
    coll_deexc_coeff = CollDeexcRateCoeff(None).calculate(
        lte_level_boltzmann_factor, coll_exc_coeff
    )
    coll_deexc_coeff.index = coll_deexc_coeff.index.swaplevel(
        "level_number_lower", "level_number_upper"
    )

    return pd.concat(
        [
            _with_ion_source_and_destination_index(coll_exc_coeff),
            _with_ion_source_and_destination_index(coll_deexc_coeff),
        ]
    )


def test_iip_cmfgen_collisional_strengths(iip_atom_data: object) -> None:
    """Verify IIP collision-strength interpolation at tabulated temperatures."""
    yg_interp, _, _ = YgInterpolator(None).calculate(
        iip_atom_data.yg_data,
        iip_atom_data.collision_data_temperatures,
        iip_atom_data.lines.index,
    )
    interpolated_yg_data = Yg(None).calculate(
        yg_interp,
        iip_atom_data.yg_data.index,
        iip_atom_data.collision_data_temperatures,
    )
    npt.assert_allclose(
        interpolated_yg_data.values,
        iip_atom_data.yg_data.values,
        atol=0,
        rtol=1e-8,
    )


def test_thermal_collision_rates_against_iip(
    iip_atom_data: object,
    iip_collision_rate_coefficients: pd.DataFrame,
) -> None:
    """Verify equilibrium collisional rates against IIP plasma properties."""
    radiative_transitions = iip_atom_data.lines.loc[
        (1, 0, slice(None), slice(None)),
    ]
    collision_strengths = iip_atom_data.yg_data.loc[
        (1, 0, slice(None), slice(None)),
    ]
    thermal_rate_solver = ThermalCollisionalRateSolver(
        iip_atom_data.levels,
        radiative_transitions,
        iip_atom_data.collision_data_temperatures,
        collision_strengths,
        collision_strengths_type="cmfgen",
        collisional_strength_approximation="regemorter",
    )

    thermal_rates = thermal_rate_solver.solve(np.array([10000.0, 20000.0]) * u.K)
    pdt.assert_frame_equal(
        thermal_rates.loc[iip_collision_rate_coefficients.index],
        iip_collision_rate_coefficients,
        check_names=False,
        check_column_type=False,
        atol=0,
        rtol=2e-5,
    )
