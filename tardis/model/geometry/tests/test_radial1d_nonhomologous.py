import numpy as np
import numpy.testing as npt
import pytest

from tardis import constants as const
from tardis.model.geometry.radial1d_nonhomologous import (
    NumbaNonhomologousRadial1DGeometry,
)

C_SPEED_OF_LIGHT = const.c.to("cm/s").value


@pytest.fixture
def numba_nonhomologous_geometry():
    return NumbaNonhomologousRadial1DGeometry(
        np.array([1.0e14, 2.0e14]),
        np.array([2.0e14, 3.0e14]),
        np.array([1.0e9, 2.0e9]),
        np.array([2.0e9, 3.0e9]),
    )


def test_numba_nonhomologous_geometry_doppler_factor(
    numba_nonhomologous_geometry,
):
    radius = 1.5e14
    mu = 0.3
    shell_id = 0
    beta = (
        numba_nonhomologous_geometry.get_velocity(radius, shell_id)
        / C_SPEED_OF_LIGHT
    )
    doppler_factor = 1.0 - mu * beta

    npt.assert_allclose(
        numba_nonhomologous_geometry.get_doppler_factor(
            radius,
            mu,
            shell_id,
            False,
        ),
        doppler_factor,
    )
    npt.assert_allclose(
        numba_nonhomologous_geometry.get_inverse_doppler_factor(
            radius,
            mu,
            shell_id,
            False,
        ),
        1.0 / doppler_factor,
    )


def test_numba_nonhomologous_geometry_doppler_factor_full_relativity(
    numba_nonhomologous_geometry,
):
    with pytest.raises(
        NotImplementedError,
        match=r"Full relativity not implemented for non-homologous mode.",
    ):
        numba_nonhomologous_geometry.get_doppler_factor(
            1.5e14,
            0.3,
            0,
            True,
        )
