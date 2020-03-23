import os
import pytest

from ctypes import (
        CDLL,
        byref,
        c_int64,
        c_double,
        c_ulong,
        )

from tardis.montecarlo import montecarlo
from tardis.montecarlo.struct import (
    RPacket, StorageModel, RKState,
    TARDIS_PACKET_STATUS_IN_PROCESS,
    CONTINUUM_OFF,
    BoundFreeTreatment
)


# Wrap the shared object containing C methods, which are tested here.
@pytest.fixture(scope='session')
def clib():
    return CDLL(os.path.join(montecarlo.__file__))


@pytest.fixture(scope="function")
def packet():
    """Fixture to return `RPacket` object with default params initialized."""
    return RPacket(
        nu=0.4,
        mu=0.3,
        energy=0.9,
        r=7.5e14,
        tau_event=2.9e13,
        nu_line=0.2,
        current_shell_id=0,
        next_line_id=1,
        last_line=0,
        close_line=0,
        current_continuum_id=1,
        virtual_packet_flag=1,
        virtual_packet=0,
        next_shell_id=1,
        status=TARDIS_PACKET_STATUS_IN_PROCESS,
        id=0,
        chi_cont=6.652486e-16,
        chi_bf_tmp_partial=(c_double * 2)(),
        compute_chi_bf=True,
        vpacket_weight=1.0
    )


@pytest.fixture(scope="function")
def model():
    """
    Fixture to return `StorageModel` object with default params initialized.
    """
    return StorageModel(
        last_line_interaction_in_id=(c_int64 * 2)(*([0] * 2)),
        last_line_interaction_shell_id=(c_int64 * 2)(*([0] * 2)),
        last_interaction_type=(c_int64 * 2)(*([2])),
        last_interaction_out_type=(c_int64 * 1)(*([0])),

        no_of_packets=2,
        no_of_shells=2,

        r_inner=(c_double * 2)(*[6.912e14, 8.64e14]),
        r_outer=(c_double * 2)(*[8.64e14, 1.0368e15]),

        time_explosion=5.2e7,
        inverse_time_explosion=1 / 5.2e7,

        electron_densities=(c_double * 2)(*[1.0e9] * 2),
        inverse_electron_densities=(c_double * 2)(*[1.0e-9] * 2),

        line_list_nu=(c_double * 5)(*[
            1.26318289e+16,
            1.26318289e+16,
            1.23357675e+16,
            1.23357675e+16,
            1.16961598e+16]),

        continuum_list_nu=(c_double * 2)(*([1.e13] * 2)),

        line_lists_tau_sobolevs=(c_double * 1000)(*([1.e-5] * 1000)),
        line_lists_j_blues=(c_double * 2)(*([1.e-10] * 2)),
        line_lists_j_blues_nd=0,

        # Init to an explicit array
        line_lists_Edotlu=(c_double * 3)(*[0.0, 0.0, 1.0]),

        no_of_lines=5,
        no_of_edges=2,

        line_interaction_id=0,
        line2macro_level_upper=(c_int64 * 2)(*([0] * 2)),

        js=(c_double * 2)(*([0.0] * 2)),
        nubars=(c_double * 2)(*([0.0] * 2)),

        spectrum_start_nu=1.e14,
        spectrum_delta_nu=293796608840.0,
        spectrum_end_nu=6.e15,

        spectrum_virt_start_nu=1e14,
        spectrum_virt_end_nu=6e15,
        spectrum_virt_nu=(c_double * 20000)(*([0.0] * 20000)),

        sigma_thomson=6.652486e-25,
        inverse_sigma_thomson=1 / 6.652486e-25,

        inner_boundary_albedo=0.0,
        reflective_inner_boundary=0,

        chi_ff_factor=(c_double * 2)(*([1.0] * 2)),
        t_electrons=(c_double * 2)(*([1.0e4] * 2)),

        l_pop=(c_double * 20000)(*([2.0] * 20000)),
        l_pop_r=(c_double * 20000)(*([3.0] * 20000)),
        cont_status=CONTINUUM_OFF,
        bf_treatment=BoundFreeTreatment.LIN_INTERPOLATION.value,
        ff_heating_estimator=(c_double * 2)(*([0.0] * 2)),
        cont_edge2macro_level=(c_int64 * 6)(*([1] * 6)),

        survival_probability=0.0,
        tau_russian=10.0,
        tau_bias=(c_double * 3)(*([5.0, 0.5, 0.0])),
        enable_biasing=0
    )

@pytest.fixture(scope="function")
def montecarlo_main_loop_model():
    """
    Fixture to return `StorageModel` object with default params initialized.
    """
    return StorageModel(
        packet_nus=(c_double * 10)(*[
            1091897396319023.375000,
            889186386678450.875000,
            1017365081318044.875000,
            1361992918139904.250000,
            429253151342138.937500,
            1107922422758092.750000,
            1486327823110125.250000,
            692412625449050.875000,
            694686330775160.625000,
            1488141901975641.000000
            ]),
        packet_mus=(c_double * 10)(*[
            0.797649,
            0.602617,
            0.510227,
            0.724160,
            0.704934,
            0.594829,
            0.647901,
            0.973652,
            0.924580,
            0.763061
            ]),
        packet_energies=(c_double * 10)(*([0.100000] * 10)),
        js=(c_double * 10)(*([0.000000] * 10)),
        nubars=(c_double * 10)(*([0.000000] * 10)),
        last_line_interaction_in_id=(c_int64 * 10)(*([-1] * 10)),
        last_line_interaction_shell_id=(c_int64 * 10)(*([-1] * 10)),
        last_interaction_type=(c_int64 * 10)(*([-1] * 10)),
        last_interaction_out_type=(c_int64 * 10)(*([5] * 10)),

        no_of_shells=20,
        no_of_packets=10,

        r_inner=(c_double * 10)(*[
            1235520000000000.000000,
            1286064000000000.000000,
            1336608000000000.000000,
            1387152000000000.000000,
            1437696000000000.000000,
            1488240000000000.000000,
            1538784000000000.000000,
            1589328000000000.000000,
            1639872000000000.000000,
            1690416000000000.000000
            ]),
        r_outer=(c_double * 10)(*[
            1286064000000000.000000,
            1336608000000000.000000,
            1387152000000000.000000,
            1437696000000000.000000,
            1488240000000000.000000,
            1538784000000000.000000,
            1589328000000000.000000,
            1639872000000000.000000,
            1690416000000000.000000,
            1740960000000000.000000
            ]),

        time_explosion=1123200.000000,
        inverse_time_explosion=1 / 1123200.000000,

        electron_densities=(c_double * 10)(*[
            2930269167.718943,
            2216660587.370028,
            1702003727.989302,
            1370776197.494959,
            1293378160.951791,
            1041706293.548131,
            641066803.873749,
            513778919.965030,
            414874363.317187,
            337506787.641974
            ]),
        inverse_electron_densities=(c_double * 10)(*[
            1 / 2930269167.718943,
            1 / 2216660587.370028,
            1 / 1702003727.989302,
            1 / 1370776197.494959,
            1 / 1293378160.951791,
            1 / 1041706293.548131,
            1 / 641066803.873749,
            1 / 513778919.965030,
            1 / 414874363.317187,
            1 / 337506787.641974
            ]),

        line_list_nu=(c_double * 10)(*[
            44897929970646376.000000,
            43320105485232064.000000,
            43300083482581312.000000,
            43167903755327728.000000,
            42116921369466568.000000,
            42062556367769000.000000,
            42045448655016688.000000,
            41909672179273904.000000,
            41893273990022496.000000,
            41126050537752416.000000
        ]),

        # continuum_list_nu=(c_double * 2)(*([1.e13] * 2)),

        line_lists_tau_sobolevs=(c_double * 10)(*([0.0] * 10)),
        line_lists_j_blues=(c_double * 10)(*([0.000000] * 10)),
        # line_lists_j_blues_nd=0,

        # # Init to an explicit array
        # line_lists_Edotlu=(c_double * 3)(*[0.0, 0.0, 1.0]),

        no_of_lines=29224,
        #no_of_edges=140309783084296,

        # line_interaction_id=0,
        # line2macro_level_upper=(c_int64 * 2)(*([0] * 2)),

        spectrum_start_nu=149896228999999.968750,
        spectrum_delta_nu=584595293100.000000,
        spectrum_end_nu=5995849159999999.000000,

        spectrum_virt_start_nu=11991698319999.998047,
        spectrum_virt_end_nu=59958491599999992.000000,
        spectrum_virt_nu=(c_double * 10)(*([0.0] * 10)),

        sigma_thomson=1 / 1503203612356297262891008.000000,
        inverse_sigma_thomson=1503203612356297262891008.000000,

        inner_boundary_albedo=0.0,
        reflective_inner_boundary=0,

        chi_ff_factor=(c_double * 10)(*([1.0] * 2)),
        t_electrons=(c_double * 10)(*([1.0e4] * 2)),

        l_pop=(c_double * 1)(*([0.0] * 1)),
        # l_pop_r=(c_double * 20000)(*([3.0] * 20000)),
        cont_status=CONTINUUM_OFF,
        bf_treatment=BoundFreeTreatment.LIN_INTERPOLATION.value,
        ff_heating_estimator=(c_double * 10)(*([0.0] * 2)),
        cont_edge2macro_level=(c_int64 * 30)(*([1] * 6)),

        survival_probability=0.0,
        tau_russian=10.000000,
        tau_bias=(c_double * 10)(*([0.0] * 10)),
        enable_biasing=0
    )

@pytest.fixture(scope="function")
def mt_state():
    """Fixture to return `RKState` object with default params initialized."""
    return RKState(
        key=(c_ulong * 624)(*([0] * 624)),
        pos=0,
        has_gauss=0,
        gauss=0.0
    )


@pytest.fixture(scope="function")
def mt_state_seeded(clib, mt_state):
    seed = 23111963
    clib.rk_seed(seed, byref(mt_state))
    return mt_state
