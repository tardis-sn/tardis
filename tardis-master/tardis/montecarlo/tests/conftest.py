import os
import pytest
from tardis.io import config_reader


from ctypes import (
    CDLL,
    byref,
    c_int64,
    c_double,
    c_ulong,
)

from tardis.montecarlo.struct import (
    RPacket,
    StorageModel,
    RKState,
    TARDIS_PACKET_STATUS_IN_PROCESS,
    CONTINUUM_OFF,
    BoundFreeTreatment,
)
from tardis.montecarlo.base import MontecarloRunner


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
        vpacket_weight=1.0,
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
        no_of_shells=2,
        r_inner=(c_double * 2)(*[6.912e14, 8.64e14]),
        r_outer=(c_double * 2)(*[8.64e14, 1.0368e15]),
        time_explosion=5.2e7,
        inverse_time_explosion=1 / 5.2e7,
        electron_densities=(c_double * 2)(*[1.0e9] * 2),
        inverse_electron_densities=(c_double * 2)(*[1.0e-9] * 2),
        line_list_nu=(c_double * 5)(
            *[
                1.26318289e16,
                1.26318289e16,
                1.23357675e16,
                1.23357675e16,
                1.16961598e16,
            ]
        ),
        continuum_list_nu=(c_double * 2)(*([1.0e13] * 2)),
        line_lists_tau_sobolevs=(c_double * 1000)(*([1.0e-5] * 1000)),
        line_lists_j_blues=(c_double * 2)(*([1.0e-10] * 2)),
        line_lists_j_blues_nd=0,
        # Init to an explicit array
        line_lists_Edotlu=(c_double * 3)(*[0.0, 0.0, 1.0]),
        no_of_lines=5,
        no_of_edges=2,
        line_interaction_id=0,
        line2macro_level_upper=(c_int64 * 2)(*([0] * 2)),
        js=(c_double * 2)(*([0.0] * 2)),
        nubars=(c_double * 2)(*([0.0] * 2)),
        spectrum_start_nu=1.0e14,
        spectrum_delta_nu=293796608840.0,
        spectrum_end_nu=6.0e15,
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
        enable_biasing=0,
    )


@pytest.fixture(scope="function")
def mt_state():
    """Fixture to return `RKState` object with default params initialized."""
    return RKState(
        key=(c_ulong * 624)(*([0] * 624)), pos=0, has_gauss=0, gauss=0.0
    )


@pytest.fixture(scope="function")
def mt_state_seeded(clib, mt_state):
    seed = 23111963
    clib.rk_seed(seed, byref(mt_state))
    return mt_state


@pytest.fixture(scope="function")
def runner():
    config_fname = "tardis/io/tests/data/tardis_configv1_verysimply.yml"
    config = config_reader.Configuration.from_yaml(config_fname)
    runner = MontecarloRunner.from_config(config)
    return runner
