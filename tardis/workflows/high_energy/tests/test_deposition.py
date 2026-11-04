from types import SimpleNamespace

import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
from astropy import units as u

from tardis.workflows.high_energy.deposition import prepare_deposition_state


def make_result() -> SimpleNamespace:
    """Return two shells and two high-energy time bins labeled in seconds."""
    return SimpleNamespace(
        gamma_ray_deposition_estimator=pd.DataFrame(
            [[10.0, 20.0], [30.0, 40.0]],
            index=["inner", "outer"],
            columns=[86400.0, 172800.0],
        ),
        positron_energy=pd.DataFrame(
            [[2.0, 8.0], [6.0, 16.0]],
            index=["inner", "outer"],
            columns=[86400.0, 172800.0],
        ),
    )


def test_deposition_selects_epoch_and_adds_positrons_once() -> None:
    """Use the containing bin and preserve the represented shell energy."""
    result = make_result()
    shell_volume = np.array([[1.0, 2.0], [3.0, 4.0]]) * u.cm**3

    deposition = prepare_deposition_state(
        result,
        np.array([1.0, 2.0, 4.0]) * u.d,
        3.0 * u.d,
        shell_volume,
        pd.Index(["inner", "outer"]),
    )

    npt.assert_array_equal(
        deposition.gamma_rate_density.to_value(u.erg / u.s / u.cm**3),
        [20.0, 40.0],
    )
    npt.assert_array_equal(
        deposition.positron_rate_density.to_value(u.erg / u.s / u.cm**3),
        [2.0 / 86400.0, 2.0 / 86400.0],
    )
    npt.assert_array_equal(
        deposition.rate_density.to_value(u.erg / u.s / u.cm**3),
        [20.0 + 2.0 / 86400.0, 40.0 + 2.0 / 86400.0],
    )
    recovered_positron_energy = (
        deposition.positron_rate_density * shell_volume[:, 1] * (2.0 * u.d)
    )
    npt.assert_array_equal(
        recovered_positron_energy.to_value(u.erg), [8.0, 16.0]
    )
    assert deposition.includes_pair_creation is False


def test_deposition_aligns_reordered_shell_labels() -> None:
    """Shell-label order does not change the physical deposition values."""
    result = make_result()
    result.gamma_ray_deposition_estimator = (
        result.gamma_ray_deposition_estimator.iloc[::-1]
    )
    result.positron_energy = result.positron_energy.iloc[::-1]

    deposition = prepare_deposition_state(
        result,
        np.array([1.0, 2.0, 4.0]) * u.d,
        1.5 * u.d,
        np.array([[1.0, 2.0], [3.0, 4.0]]) * u.cm**3,
        pd.Index(["inner", "outer"]),
    )

    npt.assert_array_equal(
        deposition.gamma_rate_density.to_value(u.erg / u.s / u.cm**3),
        [10.0, 30.0],
    )


def test_deposition_rejects_mismatched_shell_grid() -> None:
    """A nonmatching high-energy shell grid fails explicitly."""
    with pytest.raises(ValueError, match="shell labels"):
        prepare_deposition_state(
            make_result(),
            np.array([1.0, 2.0, 4.0]) * u.d,
            1.5 * u.d,
            np.ones((2, 2)) * u.cm**3,
            pd.Index(["inner", "missing"]),
        )
