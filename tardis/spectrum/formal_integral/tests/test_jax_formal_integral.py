import os

import numpy as np
import numpy.testing as ntest
import pandas as pd
import pytest

import jax
import jax.numpy as jnp

from tardis.spectrum.formal_integral.formal_integral_jax import reverse_binary_search_jax, line_search_jax, intensity_black_body_jax, populate_z_jax, calculate_z_jax
from tardis.util.base import intensity_black_body
from tardis.spectrum.formal_integral.base import C_INV
from tardis.spectrum.formal_integral.formal_integral_numba import populate_z as populate_z_numba
from tardis import constants as c
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry

TESTDATA = [
    {
        "r": np.linspace(1, 2, 3, dtype=np.float64),
    },
    {
        "r": np.linspace(0, 1, 3),
    }
]
@pytest.fixture(scope="function", params=TESTDATA)
def formal_integral_geometry(request):
    r = request.param["r"]
    geometry = NumbaRadial1DGeometry(
        r[:-1],
        r[1:],
        r[:-1] * c.c.cgs.value,
        r[1:] * c.c.cgs.value,
    )
    return geometry

@pytest.fixture
def geom(simulation_verysimple):
    # store the r_inner/p
    return simulation_verysimple.simulation_state.geometry

@pytest.fixture
def texp(simulation_verysimple):
    return simulation_verysimple.simulation_state.time_explosion.value

@pytest.fixture
def tau_sobolev(simulation_verysimple):
    return simulation_verysimple.transport.transport_state.opacity_state.tau_sobolev

@pytest.fixture
def tau_sizes(tau_sobolev):
    return tau_sobolev.shape

@pytest.fixture
def line_list_nu(simulation_verysimple):
    return simulation_verysimple.transport.transport_state.opacity_state.line_list_nu

@pytest.fixture
def mct_state(simulation_verysimple):
    return simulation_verysimple.transport.transport_state.montecarlo_state

@pytest.fixture
def electron_densities(simulation_verysimple):
    return simulation_verysimple.transport.transport_state.opacity_state.electron_densities


@pytest.mark.parametrize(
    ["x", "x_insert", "imin", "imax", "expected_params"],
    [
        (
            jnp.array([5.0, 4.0, 3.0, 1.0]),
            2.0,
            0,
            3,
            {"result": 2},
        ),
        # ( # TODO: figure out how to bound check or if should even try
        #     jnp.array([5.0, 4.0, 3.0, 2.0]),
        #     0.0,
        #     0,
        #     3,
        #     {"result": 0},  # This one might need to check for a bounds error
        # ),
    ],
)
def test_reverse_binary_search(x, x_insert, imin, imax, expected_params):
    obtained_result = 0

    obtained_result = reverse_binary_search_jax(
        x, x_insert, imin, imax
    )

    assert obtained_result == expected_params["result"]


@pytest.mark.parametrize(
    ["nu", "nu_insert", "expected_params"],
    [
        (
            jnp.array([0.5, 0.4, 0.3, 0.1]),
            0.2,
            {"result": 3},
        ),
        (
            jnp.array([0.5, 0.4, 0.3, 0.2]),
            0.1,
            {"result": 4},
        ),
        (
            jnp.array([0.4, 0.3, 0.2, 0.1]),
            0.5,
            {"result": 0},
        ),
    ],
)
def test_line_search(nu, nu_insert, expected_params):
    obtained_result = 0

    obtained_result = line_search_jax(
        nu, nu_insert
    )

    assert obtained_result == expected_params["result"]


@pytest.mark.parametrize(
    ["nu", "temperature"],
    [
        (1e14, 1e4),
        (0, 1),
        (1, 1),
    ],
)
def test_intensity_black_body(nu, temperature):
    actual = intensity_black_body_jax(nu, temperature)
    # print(actual, type(actual))
    expected = intensity_black_body(nu, temperature)
    ntest.assert_almost_equal(actual, expected)
      
@pytest.mark.parametrize("p", [0.0, 0.5, 1.0])
def test_calculate_z(formal_integral_geometry, texp, p):

    inv_t = 1.0 / texp
    r_outer = formal_integral_geometry.r_outer
    
    for r in r_outer:
        actual = calculate_z_jax(r, p, inv_t)
        if p >= r:
            assert actual == 0
        else:
            desired = np.sqrt(r * r - p * p) * C_INV * inv_t
            ntest.assert_almost_equal(actual, desired)


@pytest.mark.parametrize("ps", [
                                [0, 0.5, 1], # in photosphere
                                [1e-5, 0.5, 0.99, 1] # outside photosphere
                                ])
def test_populate_z(ps, tau_sizes, texp, formal_integral_geometry):

    _, size_shell = tau_sizes

    # populate z with the Numba formal integral
    zs_n = np.zeros((len(ps), 2 * size_shell), dtype=np.float64)
    shell_ids_n = np.zeros((len(ps), 2 * size_shell), dtype=np.int64)
    sizes_n = np.zeros(len(ps), dtype=np.int64)
    for i, p in enumerate(ps): 
        z = np.zeros(2 * size_shell, dtype=np.float64)  
        shell_id = np.zeros(2 * size_shell, dtype=np.int64)
        size_z = populate_z_numba(formal_integral_geometry, texp, p, z, shell_id)
        # print(p, z)
        zs_n[i] = z
        sizes_n[i] = size_z
        shell_ids_n[i] = shell_id

    # populate zs with the JAX formal integral
    zs_j, shell_ids_j, sizes_j = populate_z_jax(jnp.array(ps), jnp.array(formal_integral_geometry.r_inner), jnp.array(formal_integral_geometry.r_outer), texp, size_shell)
    
    assert np.allclose(zs_n, zs_j)
    assert np.allclose(shell_ids_n, shell_ids_j)
    assert np.allclose(sizes_n, sizes_j)