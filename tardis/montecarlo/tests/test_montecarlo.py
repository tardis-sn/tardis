"""
Unit tests for methods in `tardis/montecarlo/src/cmontecarlo.c`.
* `ctypes` library is used to wrap C methods and expose them to python.


Probable Reasons for Failing Tests:
-----------------------------------

1. Change made in C struct declarations:
  - Reflect the changes done in C structs, into Python counterparts.
  - Check **tardis/montecarlo/struct.py**.

2. Return type of any method changed:
  - Modify the `restype` parameter in the test method here.
  - For example:
        ```
        clib.rpacket_doppler_factor.restype = c_double
        ```

3. Underlying logic modified:
  - Check whether the changes made in C method are logically correct.
  - If the changes made were correct and necessary, update the corresponding
    test case.


General Test Design Procedure:
------------------------------

Please follow this design procedure while adding a new test:

1. Parametrization as per desire of code coverage.
  - C tests have different flows controlled by conditional statements.
    Parameters checked in conditions can be provided in different testcases.
  - Keep consistency with variable names as (in order):
    - `packet_params`
    - `model_params`
    - `expected_params` (`expected` if only one value to be asserted.)
  - Suggested variable names can be compromised if readability of the test
    increases.

2. Test Method body:
  - Keep name as `test_` + `(name of C method)`.
  - Refer to method `test_rpacket_doppler_factor` below for description.
"""


import os
import pytest
import numpy as np
import pandas as pd
import tardis.montecarlo.formal_integral as formal_integral
import tardis.montecarlo.montecarlo_numba.r_packet as r_packet
import tardis.montecarlo.montecarlo_configuration as mc
from tardis import constants as const
from tardis.montecarlo.montecarlo_numba.numba_interface import Estimators
from tardis.montecarlo.montecarlo_numba import macro_atom

C_SPEED_OF_LIGHT = const.c.to("cm/s").value

pytestmark = pytest.mark.skip(reason="Port from C to numba")


from ctypes import (
    CDLL,
    byref,
    c_uint,
    c_int64,
    c_double,
    c_ulong,
    c_void_p,
    cast,
    POINTER,
    pointer,
    CFUNCTYPE,
)

from numpy.testing import (
    assert_equal,
    assert_almost_equal,
    assert_array_equal,
    assert_allclose,
)

from tardis import __path__ as path
from tardis.montecarlo.struct import (
    RPacket,
    StorageModel,
    RKState,
    PhotoXsect1level,
    TARDIS_ERROR_OK,
    TARDIS_ERROR_BOUNDS_ERROR,
    TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE,
    TARDIS_PACKET_STATUS_IN_PROCESS,
    TARDIS_PACKET_STATUS_EMITTED,
    TARDIS_PACKET_STATUS_REABSORBED,
    CONTINUUM_OFF,
    CONTINUUM_ON,
    INVERSE_C,
    BoundFreeTreatment,
)


@pytest.fixture(scope="module")
def continuum_compare_data_fname():
    fname = "continuum_compare_data.hdf"
    return os.path.join(path[0], "montecarlo", "tests", "data", fname)


@pytest.fixture(scope="module")
def continuum_compare_data(continuum_compare_data_fname, request):
    compare_data = pd.HDFStore(continuum_compare_data_fname, mode="r")

    def fin():
        compare_data.close()

    request.addfinalizer(fin)

    return compare_data


@pytest.fixture(scope="function")
def expected_ff_emissivity(continuum_compare_data):
    emissivities = continuum_compare_data["ff_emissivity"]

    def ff_emissivity(t_electron):
        emissivity = emissivities[t_electron]
        nu_bins = emissivity["nu_bins"].values
        emissivity_value = emissivity["emissivity"].dropna().values

        return nu_bins, emissivity_value

    return ff_emissivity


@pytest.fixture(scope="module")
def get_rkstate(continuum_compare_data):
    data = continuum_compare_data["z2rkstate_key"]
    pos_data = continuum_compare_data["z2rkstate_pos"]

    def z2rkstate(z_random):
        key = (c_ulong * 624)(*data.loc[z_random].values)
        pos = pos_data.loc[z_random]
        return RKState(key=key, pos=pos, has_gauss=0, gauss=0.0)

    return z2rkstate


@pytest.fixture(scope="function")
def model_w_edges(ion_edges, model):
    photo_xsect = (POINTER(PhotoXsect1level) * len(ion_edges))()

    for i, edge in enumerate(ion_edges):
        x_sect_1level = PhotoXsect1level()
        for key, value in edge.items():
            if key in ["nu", "x_sect"]:
                value = (c_double * len(value))(*value)
            setattr(x_sect_1level, key, value)
        photo_xsect[i] = pointer(x_sect_1level)

    no_of_edges = len(ion_edges)
    continuum_list_nu = (c_double * no_of_edges)(
        *[edge["nu"][0] for edge in ion_edges]
    )

    model.photo_xsect = photo_xsect
    model.continuum_list_nu = continuum_list_nu
    model.no_of_edges = no_of_edges

    estimator_size = model.no_of_shells * no_of_edges
    estims = [
        "photo_ion_estimator",
        "stim_recomb_estimator",
        "bf_heating_estimator",
        "stim_recomb_cooling_estimator",
    ]
    for estimator in estims:
        setattr(
            model, estimator, (c_double * estimator_size)(*[0] * estimator_size)
        )

    model.photo_ion_estimator_statistics = (c_int64 * estimator_size)(
        *[0] * estimator_size
    )
    return model


@pytest.fixture(scope="module")
def ion_edges():
    return [
        {
            "nu": [4.0e14, 4.1e14, 4.2e14, 4.3e14],
            "x_sect": [1.0, 0.9, 0.8, 0.7],
            "no_of_points": 4,
        },
        {
            "nu": [3.0e14, 3.1e14, 3.2e14, 3.3e14, 3.4e14],
            "x_sect": [1.0, 0.9, 0.8, 0.7, 0.6],
            "no_of_points": 5,
        },
        {
            "nu": [2.8e14, 3.0e14, 3.2e14, 3.4e14],
            "x_sect": [2.0, 1.8, 1.6, 1.4],
            "no_of_points": 4,
        },
    ]


@pytest.fixture(scope="module")
def mock_sample_nu():
    SAMPLE_NUFUNC = CFUNCTYPE(
        c_double, POINTER(RPacket), POINTER(StorageModel), POINTER(RKState)
    )

    def sample_nu_simple(packet, model, mt_state):
        return packet.contents.nu

    return SAMPLE_NUFUNC(sample_nu_simple)


@pytest.fixture(scope="function")
def model_3lvlatom(model):
    model.line2macro_level_upper = (c_int64 * 3)(*[2, 1, 2])
    model.macro_block_references = (c_int64 * 3)(*[0, 2, 5])

    transition_probabilities = [
        0.0,
        0.0,
        0.75,
        0.25,
        0.0,
        0.25,
        0.5,
        0.25,
        0.0,  # shell_id = 0
        0.0,
        0.0,
        1.00,
        0.00,
        0.0,
        0.00,
        0.0,
        1.00,
        0.0,  # shell_id = 1
    ]

    nd = len(transition_probabilities) // 2
    model.transition_type = (c_int64 * nd)(*[1, 1, -1, 1, 0, 0, -1, -1, 0])
    model.destination_level_id = (c_int64 * nd)(*[1, 2, 0, 2, 0, 1, 1, 0, 0])
    model.transition_line_id = (c_int64 * nd)(*[0, 1, 1, 2, 1, 2, 2, 0, 0])

    model.transition_probabilities_nd = c_int64(nd)
    model.transition_probabilities = (c_double * (nd * 2))(
        *transition_probabilities
    )

    model.last_line_interaction_out_id = (c_int64 * 1)(*[-5])

    return model


def d_cont_setter(d_cont, model, packet):
    model.inverse_electron_densities[packet.current_shell_id] = c_double(1.0)
    model.inverse_sigma_thomson = c_double(1.0)
    packet.tau_event = c_double(d_cont)


def d_line_setter(d_line, model, packet):
    packet.mu = c_double(0.0)
    scale = d_line * 1e1
    model.time_explosion = c_double(INVERSE_C * scale)
    packet.nu = c_double(1.0)
    nu_line = 1.0 - d_line / scale
    packet.nu_line = c_double(nu_line)


def d_boundary_setter(d_boundary, model, packet):
    packet.mu = c_double(1e-16)
    r_outer = 2.0 * d_boundary
    model.r_outer[packet.current_shell_id] = r_outer

    r = np.sqrt(r_outer ** 2 - d_boundary ** 2)
    packet.r = r


"""
Important Tests:
----------------
The tests written further (till next block comment is encountered) have been
categorized as important tests, these tests correspond to methods which are
relatively old and stable code.
"""


@pytest.mark.parametrize(
    ["x", "x_insert", "imin", "imax", "expected_params"],
    [
        (
            [5.0, 4.0, 3.0, 1.0],
            2.0,
            0,
            3,
            {"result": 2, "ret_val": TARDIS_ERROR_OK},
        ),
        (
            [5.0, 4.0, 3.0, 2.0],
            0.0,
            0,
            3,
            {"result": 0, "ret_val": TARDIS_ERROR_BOUNDS_ERROR},
        ),
    ],
)
def test_reverse_binary_search(x, x_insert, imin, imax, expected_params):
    # x = (c_double * (imax - imin + 1))(*x)
    obtained_result = 0

    obtained_tardis_error = formal_integral.reverse_binary_search(
        x, x_insert, imin, imax, obtained_result
    )

    assert obtained_result.value == expected_params["result"]
    assert obtained_tardis_error == expected_params["ret_val"]


@pytest.mark.parametrize(
    ["nu", "nu_insert", "number_of_lines", "expected_params"],
    [
        (
            [0.5, 0.4, 0.3, 0.1],
            0.2,
            4,
            {"result": 3, "ret_val": TARDIS_ERROR_OK},
        ),
        (
            [0.5, 0.4, 0.3, 0.2],
            0.1,
            4,
            {"result": 4, "ret_val": TARDIS_ERROR_OK},
        ),
        (
            [0.4, 0.3, 0.2, 0.1],
            0.5,
            4,
            {"result": 0, "ret_val": TARDIS_ERROR_OK},
        ),
    ],
)
def test_line_search(nu, nu_insert, number_of_lines, expected_params):
    # nu = (c_double * number_of_lines)(*nu)
    obtained_result = 0

    obtained_tardis_error = formal_integral.line_search(
        nu, nu_insert, number_of_lines, obtained_result
    )

    assert obtained_result.value == expected_params["result"]
    assert obtained_tardis_error == expected_params["ret_val"]


@pytest.mark.parametrize(
    ["x", "x_insert", "imin", "imax", "expected_params"],
    [
        (
            [2.0, 4.0, 6.0, 7.0],
            5.0,
            0,
            3,
            {"result": 2, "ret_val": TARDIS_ERROR_OK},
        ),
        (
            [2.0, 3.0, 5.0, 7.0],
            8.0,
            0,
            3,
            {"result": 0, "ret_val": TARDIS_ERROR_BOUNDS_ERROR},
        ),
        (
            [2.0, 4.0, 6.0, 7.0],
            4.0,
            0,
            3,
            {"result": 0, "ret_val": TARDIS_ERROR_OK},
        ),
    ],
)
def test_binary_search(x, x_insert, imin, imax, expected_params):
    obtained_result = 0

    obtained_tardis_error = formal_integral.binary_search(
        x, x_insert, imin, imax, obtained_result
    )

    assert obtained_result.value == expected_params["result"]
    assert obtained_tardis_error == expected_params["ret_val"]


@pytest.mark.parametrize(
    ["mu", "r", "inv_t_exp", "expected"],
    [
        (0.3, 7.5e14, 1 / 5.2e7, 0.9998556693818854),
        (-0.3, 0, 1 / 2.6e7, 1.0),
        (0, 1, 1 / 2.6e7, 1.0),
    ],
)
def test_get_doppler_factor(mu, r, inv_t_exp, expected, packet, model):
    # Set the params from test cases here
    # TODO: add relativity tests
    time_explosion = 1 / inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables

    obtained = r_packet.get_doppler_factor(r, mu, time_explosion)

    # Perform required assertions
    assert_almost_equal(obtained, expected)


@pytest.mark.parametrize(["mu", "r", "inv_t_exp"], [(-0.3, 5, 1e10)])
def test_unphysical_doppler_factor(mu, r, inv_t_exp, packet, model):
    # Set the params from test cases here
    # TODO: add relativity tests
    time_explosion = 1 / inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables
    with pytest.raises(r_packet.SuperluminalError):
        obtained = r_packet.get_doppler_factor(r, mu, time_explosion)


@pytest.mark.parametrize(
    ["mu", "r", "inv_t_exp", "expected"],
    [
        (0.3, 7.5e14, 1 / 5.2e7, 1 / 0.9998556693818854),
        (-0.3, 0, 1 / 2.6e7, 1.0),
        (0, 1, 1 / 2.6e7, 1.0),
    ],
)
def test_get_inverse_doppler_factor(mu, r, inv_t_exp, expected, packet, model):
    # Set the params from test cases here
    # TODO: add relativity tests
    time_explosion = 1 / inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables

    obtained = r_packet.get_inverse_doppler_factor(r, mu, time_explosion)

    # Perform required assertions
    assert_almost_equal(obtained, expected)


@pytest.mark.parametrize(["mu", "r", "inv_t_exp"], [(-0.3, 5, 1e10)])
def test_unphysical_inverse_doppler_factor(mu, r, inv_t_exp, packet, model):
    # Set the params from test cases here
    # TODO: add relativity tests
    time_explosion = 1 / inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables
    with pytest.raises(r_packet.SuperluminalError):
        obtained = r_packet.get_inverse_doppler_factor(r, mu, time_explosion)


@pytest.mark.parametrize(
    ["mu", "r", "inv_t_exp", "expected"],
    [
        (-0.3, 10000000, 0.001, 1.0000001000692842),
        (-0.3, 0, 1 / 2.6e7, 1.0),
        (0, 1, 1 / 2.6e7, 1.0),
    ],
)
def test_get_doppler_factor_full_relativity(
    mu, r, inv_t_exp, expected, packet, model
):
    # Set the params from test cases here
    # TODO: add relativity tests
    mc.full_relativity = True
    time_explosion = 1 / inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables

    obtained = r_packet.get_doppler_factor(r, mu, time_explosion)
    mc.full_relativity = False
    # Perform required assertions
    assert_almost_equal(obtained, expected)


@pytest.mark.parametrize(
    ["mu", "r", "inv_t_exp", "expected"],
    [
        (-0.3, 10000000, 0.001, 0.999999899930827),
        (-0.3, 0, 1 / 2.6e7, 1.0),
        (0, 1, 1 / 2.6e7, 1.0),
    ],
)
def test_get_inverse_doppler_factor_full_relativity(
    mu, r, inv_t_exp, expected, packet, model
):
    # Set the params from test cases here
    # TODO: add relativity tests
    mc.full_relativity = True
    time_explosion = 1 / inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables

    obtained = r_packet.get_inverse_doppler_factor(r, mu, time_explosion)
    mc.full_relativity = False
    # Perform required assertions
    assert_almost_equal(obtained, expected)


def test_get_random_mu_different_output():
    """
    Ensure that different calls results
    """
    output1 = r_packet.get_random_mu()
    output2 = r_packet.get_random_mu()
    assert output1 != output2


def test_get_random_mu_different_output():
    """
    Ensure that different calls results
    """
    output1 = r_packet.get_random_mu()
    output2 = r_packet.get_random_mu()
    assert output1 != output2


@pytest.mark.parametrize(
    ["mu", "r", "time_explosion"], [(1, C_SPEED_OF_LIGHT, 1)]
)
def test_angle_ab_LF_to_CMF_diverge(mu, r, time_explosion, packet):
    """
    This equation should diverage as costheta --> 1 and beta --> 1
    """
    packet = r_packet.RPacket(packet.r, packet.mu, packet.nu, packet.energy)
    packet.r = r
    packet.mu = mu
    with pytest.raises(ZeroDivisionError):
        obtained = r_packet.angle_aberration_LF_to_CMF(
            packet, time_explosion, mu
        )


@pytest.mark.parametrize(
    ["mu", "r", "time_explosion"], [(0.3, 1e7, 1e10), (-0.3, 1e7, 1e11)]
)
def test_both_angle_aberrations(mu, r, time_explosion, packet):
    """
    The angle aberration functions should be the functional inverse of one
    another.
    """
    packet = r_packet.RPacket(r, packet.mu, packet.nu, packet.energy)
    packet.r = r
    obtained_mu = r_packet.angle_aberration_LF_to_CMF(
        packet, time_explosion, mu
    )
    inverse_obtained_mu = r_packet.angle_aberration_CMF_to_LF(
        packet, time_explosion, obtained_mu
    )
    assert_almost_equal(inverse_obtained_mu, mu)


@pytest.mark.parametrize(
    ["mu", "r", "time_explosion"], [(0.3, 7.5e14, 5.2e5), (-0.3, 7.5e14, 5.2e5)]
)
def test_both_angle_aberrations_inverse(mu, r, time_explosion, packet):
    """
    The angle aberration functions should be the functional inverse of one
    another.
    """
    packet = r_packet.RPacket(r, packet.mu, packet.nu, packet.energy)
    packet.r = r
    obtained_mu = r_packet.angle_aberration_CMF_to_LF(
        packet, time_explosion, mu
    )
    inverse_obtained_mu = r_packet.angle_aberration_LF_to_CMF(
        packet, time_explosion, obtained_mu
    )
    assert_almost_equal(inverse_obtained_mu, mu)


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, 11, 132), (132, 1, 133), (132, 2, 133)],
)
def test_move_packet_across_shell_boundary_emitted(
    packet, current_shell_id, delta_shell, no_of_shells
):
    packet = r_packet.RPacket(packet.r, packet.mu, packet.nu, packet.energy)
    packet.current_shell_id = current_shell_id
    r_packet.move_packet_across_shell_boundary(
        packet, delta_shell, no_of_shells
    )
    assert packet.status == r_packet.PacketStatus.EMITTED


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, -133, 132), (132, -133, 133), (132, -1e9, 133)],
)
def test_move_packet_across_shell_boundary_reabsorbed(
    packet, current_shell_id, delta_shell, no_of_shells
):
    packet = r_packet.RPacket(packet.r, packet.mu, packet.nu, packet.energy)
    packet.current_shell_id = current_shell_id
    r_packet.move_packet_across_shell_boundary(
        packet, delta_shell, no_of_shells
    )
    assert packet.status == r_packet.PacketStatus.REABSORBED


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, -1, 199), (132, 0, 132), (132, 20, 154)],
)
def test_move_packet_across_shell_boundary_increment(
    packet, current_shell_id, delta_shell, no_of_shells
):
    packet = r_packet.RPacket(packet.r, packet.mu, packet.nu, packet.energy)
    packet.current_shell_id = current_shell_id
    r_packet.move_packet_across_shell_boundary(
        packet, delta_shell, no_of_shells
    )
    assert packet.current_shell_id == current_shell_id + delta_shell


@pytest.mark.parametrize(
    ["distance_trace", "time_explosion", "mu", "r"],
    [(0, 1, 0, 0), (0, 1, 1, 0), (0, 1, 0, 1)],
)
def test_packet_energy_limit_one(packet, distance_trace, time_explosion, mu, r):
    initial_energy = packet.energy
    packet = r_packet.RPacket(r, mu, packet.nu, packet.energy)
    new_energy = r_packet.calc_packet_energy(
        packet, distance_trace, time_explosion
    )
    assert new_energy == initial_energy


@pytest.mark.parametrize(
    ["packet_params", "expected_params"],
    [
        ({"mu": 0.3, "r": 7.5e14}, {"d_boundary": 259376919351035.88}),
        ({"mu": -0.3, "r": 7.5e13}, {"d_boundary": -664987228972291.5}),
        ({"mu": -0.3, "r": 7.5e14}, {"d_boundary": 709376919351035.9}),
    ],
)
def test_compute_distance2boundary(
    packet_params, expected_params, packet, model
):
    mu = packet_params["mu"]
    r = packet_params["r"]

    d_boundary = r_packet.calculate_distance_boundary(
        r, mu, model.r_inner[0], model.r_outer[0]
    )

    assert_almost_equal(d_boundary[0], expected_params["d_boundary"])


#
#
# TODO: split this into two tests - one to assert errors and other for d_line
@pytest.mark.parametrize(
    ["packet_params", "expected_params"],
    [
        (
            {"nu_line": 0.1, "next_line_id": 0, "last_line": 1},
            {"tardis_error": None, "d_line": 1e99},
        ),
        (
            {"nu_line": 0.2, "next_line_id": 1, "last_line": 0},
            {"tardis_error": None, "d_line": 7.792353908000001e17},
        ),
        (
            {"nu_line": 0.5, "next_line_id": 1, "last_line": 0},
            {"tardis_error": r_packet.MonteCarloException, "d_line": 0.0},
        ),
        (
            {"nu_line": 0.6, "next_line_id": 0, "last_line": 0},
            {"tardis_error": r_packet.MonteCarloException, "d_line": 0.0},
        ),
    ],
)
def test_compute_distance2line(packet_params, expected_params, packet, model):
    packet = r_packet.RPacket(packet.r, packet.mu, packet.nu, packet.energy)
    nu_line = packet_params["nu_line"]
    # packet.next_line_id = packet_params['next_line_id']
    # packet.last_line = packet_params['last_line']

    time_explosion = model.time_explosion

    doppler_factor = r_packet.get_doppler_factor(
        packet.r, packet.mu, time_explosion
    )
    comov_nu = packet.nu * doppler_factor

    d_line = 0
    obtained_tardis_error = None
    try:
        d_line = r_packet.calculate_distance_line(
            packet, comov_nu, nu_line, time_explosion
        )
    except r_packet.MonteCarloException:
        obtained_tardis_error = r_packet.MonteCarloException

    assert_almost_equal(d_line, expected_params["d_line"])
    assert obtained_tardis_error == expected_params["tardis_error"]


# @pytest.mark.parametrize(
#     ['packet_params', 'expected_params'],
#     [({'virtual_packet': 0},
#      {'chi_cont': 6.652486e-16, 'd_cont': 4.359272608766106e+28}),
#
#      ({'virtual_packet': 1},
#       {'chi_cont': 6.652486e-16, 'd_cont': 1e+99})]
# )
# def test_compute_distance2continuum(clib, packet_params, expected_params, packet, model):
#     packet.virtual_packet = packet_params['virtual_packet']
#
#     clib.compute_distance2continuum(byref(packet), byref(model))
#
#     assert_almost_equal(packet.chi_cont, expected_params['chi_cont'])
#     assert_almost_equal(packet.d_cont, expected_params['d_cont'])


@pytest.mark.parametrize("full_relativity", [1, 0])
@pytest.mark.parametrize(
    ["packet_params", "expected_params"],
    [
        (
            {"nu": 0.4, "mu": 0.3, "energy": 0.9, "r": 7.5e14},
            {
                "mu": 0.3120599529139568,
                "r": 753060422542573.9,
                "j": 8998701024436.969,
                "nubar": 3598960894542.354,
            },
        ),
        (
            {"nu": 0.6, "mu": -0.5, "energy": 0.5, "r": 8.1e14},
            {
                "mu": -0.4906548373534084,
                "r": 805046582503149.2,
                "j": 5001298975563.031,
                "nubar": 3001558973156.1387,
            },
        ),
    ],
)
def test_move_packet(
    packet_params, expected_params, packet, model, full_relativity
):
    distance = 1e13
    packet.nu = packet_params["nu"]
    packet.mu = packet_params["mu"]
    packet.energy = packet_params["energy"]
    packet.r = packet_params["r"]
    # model.full_relativity = full_relativity
    mc.full_relativity = full_relativity

    doppler_factor = r_packet.get_doppler_factor(
        packet.r, packet.mu, model.time_explosion
    )
    numba_estimator = Estimators(
        packet_params["j"], packet_params["nu_bar"], 0, 0
    )
    r_packet.move_r_packet(
        packet, distance, model.time_explosion, numba_estimator
    )

    assert_almost_equal(packet.mu, expected_params["mu"])
    assert_almost_equal(packet.r, expected_params["r"])

    expected_j = expected_params["j"]
    expected_nubar = expected_params["nubar"]
    if full_relativity:
        expected_j *= doppler_factor
        expected_nubar *= doppler_factor

    mc.full_relativity = False

    assert_allclose(
        numba_estimator.j_estimator[packet.current_shell_id],
        expected_j,
        rtol=5e-7,
    )
    assert_allclose(
        numba_estimator.nu_bar_estimator[packet.current_shell_id],
        expected_nubar,
        rtol=5e-7,
    )


# @pytest.mark.continuumtest
# @pytest.mark.parametrize(
#     ['packet_params', 'j_blue_idx', 'expected'],
#     [({'nu': 0.30, 'energy': 0.30}, 0, 1.0),
#      ({'nu': 0.20, 'energy': 1.e5}, 0, 5e5),
#      ({'nu': 2e15, 'energy': 0.50}, 1, 2.5e-16),
#      ({'nu': 0.40, 'energy': 1e-7}, 1, 2.5e-7)],
# )
# def test_increment_j_blue_estimator_full_relativity(packet_params,
#                                                     j_blue_idx, expected,
#                                                     packet, model):
#     packet.nu = packet_params['nu']
#     packet.energy = packet_params['energy']
#     model.full_relativity = True
#
#     r_packet.increment_j_blue_estimator(byref(packet), byref(model),
#                                     c_double(packet.d_line),
#                                     c_int64(j_blue_idx))
#
#     assert_almost_equal(model.line_lists_j_blues[j_blue_idx], expected)
#
#
# @pytest.mark.parametrize(
#     ['packet_params', 'cur_line_id', 'expected'],
#     [({'nu': 0.1, 'mu': 0.3, 'r': 7.5e14}, 0, 8.998643292289723),
#      ({'nu': 0.2, 'mu': -.3, 'r': 7.7e14}, 0, 4.499971133976377),
#      ({'nu': 0.5, 'mu': 0.5, 'r': 7.9e14}, 1, 0.719988453650551),
#      ({'nu': 0.6, 'mu': -.5, 'r': 8.1e14}, 1, 0.499990378058792)]
# )
# def test_increment_j_blue_estimator(packet_params, cur_line_id, expected, packet):
#
#     numba_interface
#     packet = r_packet.RPacket(packet_params['r'],
#                               packet_params['mu'],
#                               packet_params['nu'],
#                               packet.energy)
#
#     r_packet.compute_distance2line(byref(packet), byref(model))
#     r_packet.move_r_packet(packet,
#                            distance,
#                            model.time_explosion,
#                            numba_estimator)
#     r_packet.move_packet(byref(packet), byref(model), c_double(1.e13))
#     r_packet.update_line_estimators(estimators, r_packet,
#                                     cur_line_id,
#                                     model.distance_trace,
#                                     model.time_explosion)
#
#     assert_almost_equal(model.line_lists_j_blues[j_blue_idx], expected)


@pytest.mark.parametrize(
    ["packet_params", "expected_params"],
    [
        (
            {"nu": 0.4, "mu": 0.3, "energy": 0.9, "r": 7.5e14},
            {"nu": 0.39974659819356556, "energy": 0.8994298459355226},
        ),
        (
            {"nu": 0.6, "mu": -0.5, "energy": 0.5, "r": 8.1e14},
            {"nu": 0.5998422620533325, "energy": 0.4998685517111104},
        ),
    ],
)
def test_montecarlo_thomson_scatter(
    clib, packet_params, expected_params, packet, model, mt_state
):
    packet.nu = packet_params["nu"]
    packet.mu = packet_params["mu"]
    packet.energy = packet_params["energy"]
    packet.r = packet_params["r"]

    clib.montecarlo_thomson_scatter(
        byref(packet), byref(model), c_double(1.0e13), byref(mt_state)
    )

    assert_almost_equal(packet.nu, expected_params["nu"])
    assert_almost_equal(packet.energy, expected_params["energy"])


@pytest.mark.parametrize(
    ["packet_params", "expected_params"],
    # TODO: Add scientifically sound test cases.
    [
        (
            {"virtual_packet": 1, "tau_event": 2.9e13, "last_line": 0},
            {"tau_event": 2.9e13, "next_line_id": 2},
        ),
        (
            {"virtual_packet": 0, "tau_event": 2.9e13, "last_line": 0},
            {"tau_event": 2.9e13, "next_line_id": 2},
        ),
        (
            {"virtual_packet": 0, "tau_event": 2.9e13, "last_line": 0},
            {"tau_event": 2.9e13, "next_line_id": 2},
        ),
    ],
)
def test_montecarlo_line_scatter(
    clib, packet_params, expected_params, packet, model, mt_state
):
    packet.virtual_packet = packet_params["virtual_packet"]
    packet.tau_event = packet_params["tau_event"]
    packet.last_line = packet_params["last_line"]

    clib.montecarlo_line_scatter(
        byref(packet), byref(model), c_double(1.0e13), byref(mt_state)
    )

    assert_almost_equal(packet.tau_event, expected_params["tau_event"])
    assert_almost_equal(packet.next_line_id, expected_params["next_line_id"])


@pytest.mark.parametrize(
    ["z_random", "packet_params", "expected"],
    [
        (
            0.22443743797312765,
            {"activation_level": 1, "shell_id": 0},
            1,
        ),  # Direct deactivation
        (
            0.78961460371187597,  # next z_random = 0.818455414618
            {"activation_level": 1, "shell_id": 0},
            0,
        ),  # Upwards jump, then deactivation
        (
            0.22443743797312765,  # next z_random = 0.545678896748
            {"activation_level": 2, "shell_id": 0},
            1,
        ),  # Downwards jump, then deactivation
        (
            0.765958602560605,  # next z_random = 0.145914243888, 0.712382380384
            {"activation_level": 1, "shell_id": 0},
            1,
        ),  # Upwards jump, downwards jump, then deactivation
        (0.22443743797312765, {"activation_level": 2, "shell_id": 1}, 0),
    ],  # Direct deactivation
)
def test_macro_atom(packet, z_random, packet_params, get_rkstate, expected):
    packet.macro_atom_activation_level = packet_params["activation_level"]
    packet = r_packet.RPacket(
        r=packet.r, mu=packet.mu, nu=packet.nu, energy=packet.energy
    )
    packet.current_shell_id = packet_params["shell_id"]
    rkstate = get_rkstate(z_random)

    macro_atom.macro_atom(packet, numba_plasma)
    obtained_line_id = model_3lvlatom.last_line_interaction_out_id[packet.id]

    assert_equal(obtained_line_id, expected)


"""
Simple Tests:
----------------
These test check very simple pieces of code still work.
"""


@pytest.mark.parametrize(
    ["packet_params", "line_idx", "expected"],
    [
        ({"energy": 0.0}, 0, 0),
        ({"energy": 1.0}, 1, 1),
        ({"energy": 0.5}, 2, 1.5),
    ],
)
@pytest.mark.skipif(True, reason="Needs rewrite to be relevant.")
def test_increment_Edotlu_estimator(
    clib, packet_params, line_idx, expected, packet, model
):
    packet.energy = packet_params["energy"]

    clib.increment_Edotlu_estimator(
        byref(packet), byref(model), c_int64(line_idx)
    )

    assert_almost_equal(model.line_lists_Edotlu[line_idx], expected)


"""
Continuum Tests:
----------------
The tests written further (till next block comment is encountered) are for the
methods related to continuum interactions.
"""


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    ["mu", "r", "inv_t_exp", "full_relativity"],
    [
        (0.8, 7.5e14, 1 / 5.2e5, 1),
        (-0.7, 7.5e14, 1 / 5.2e5, 1),
        (0.3, 7.5e14, 1 / 2.2e5, 1),
        (0.0, 7.5e14, 1 / 2.2e5, 1),
        (-0.7, 7.5e14, 1 / 5.2e5, 0),
    ],
)
def test_frame_transformations(packet, mu, r, inv_t_exp, full_relativity):
    packet = r_packet.RPacket(r=r, mu=mu, energy=packet.energy, nu=packet.nu)
    mc.full_relativity = bool(full_relativity)
    mc.full_relativity = full_relativity

    inverse_doppler_factor = r_packet.get_inverse_doppler_factor(
        r, mu, 1 / inv_t_exp
    )
    r_packet.angle_aberration_CMF_to_LF(packet, 1 / inv_t_exp, packet.mu)

    doppler_factor = r_packet.get_doppler_factor(r, mu, 1 / inv_t_exp)
    mc.full_relativity = False

    assert_almost_equal(doppler_factor * inverse_doppler_factor, 1.0)


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    ["mu", "r", "inv_t_exp"],
    [
        (0.8, 7.5e14, 1 / 5.2e5),
        (-0.7, 7.5e14, 1 / 5.2e5),
        (0.3, 7.5e14, 1 / 2.2e5),
        (0.0, 7.5e14, 1 / 2.2e5),
        (-0.7, 7.5e14, 1 / 5.2e5),
    ],
)
def test_angle_transformation_invariance(packet, model, mu, r, inv_t_exp):
    packet = r_packet.RPacket(r, mu, packet.nu, packet.energy)
    model.full_relativity = 1

    mu1 = r_packet.angle_aberration_CMF_to_LF(packet, 1 / inv_t_exp, mu)
    mu_obtained = r_packet.angle_aberration_LF_to_CMF(
        packet, 1 / inv_t_exp, mu1
    )

    assert_almost_equal(mu_obtained, mu)


@pytest.mark.continuumtest
@pytest.mark.parametrize("full_relativity", [1, 0])
@pytest.mark.parametrize(
    ["mu", "r", "t_exp", "nu", "nu_line"],
    [
        (0.8, 7.5e14, 5.2e5, 1.0e15, 9.4e14),
        (0.0, 6.3e14, 2.2e5, 6.0e12, 5.8e12),
        (1.0, 9.0e14, 2.2e5, 4.0e8, 3.4e8),
        (0.9, 9.0e14, 0.5e5, 1.0e15, 4.5e14),
        (-0.7, 7.5e14, 5.2e5, 1.0e15, 9.8e14),
        (-1.0, 6.3e14, 2.2e5, 6.0e12, 6.55e12),
    ],
)
def test_compute_distance2line_relativistic(
    mu, r, t_exp, nu, nu_line, full_relativity, packet, runner
):
    packet = r_packet.RPacket(r=r, nu=nu, mu=mu, energy=packet.energy)
    # packet.nu_line = nu_line
    numba_estimator = Estimators(
        runner.j_estimator,
        runner.nu_bar_estimator,
        runner.j_blue_estimator,
        runner.Edotlu_estimator,
    )
    mc.full_relativity = bool(full_relativity)

    doppler_factor = r_packet.get_doppler_factor(r, mu, t_exp)
    comov_nu = packet.nu * doppler_factor
    distance = r_packet.calculate_distance_line(
        packet, comov_nu, nu_line, t_exp
    )
    r_packet.move_r_packet(packet, distance, t_exp, numba_estimator)

    doppler_factor = r_packet.get_doppler_factor(r, mu, t_exp)
    comov_nu = packet.nu * doppler_factor
    mc.full_relativity = False

    assert_allclose(comov_nu, nu_line, rtol=1e-14)
