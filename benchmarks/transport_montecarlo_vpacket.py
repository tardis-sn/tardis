"""
Basic TARDIS Benchmark.
"""

import numpy as np

import tardis.transport.montecarlo.vpacket as vpacket
from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.frame_transformations import (
    get_doppler_factor,
)


class BenchmarkMontecarloMontecarloNumbaVpacket(BenchmarkBase):
    """
    Class to benchmark the single packet loop function.
    """

    @property
    def v_packet(self):
        return vpacket.VPacket(
            r=7.5e14,
            nu=4e15,
            mu=0.3,
            energy=0.9,
            current_shell_id=0,
            next_line_id=0,
            index=0,
        )

    def v_packet_initialize_line_id(
        self, v_packet, opacity_state, time_explosion, enable_full_relativity
    ):
        inverse_line_list_nu = opacity_state.line_list_nu[::-1]
        doppler_factor = get_doppler_factor(
            v_packet.r, v_packet.mu, time_explosion, enable_full_relativity
        )
        comov_nu = v_packet.nu * doppler_factor
        next_line_id = len(opacity_state.line_list_nu) - np.searchsorted(
            inverse_line_list_nu, comov_nu
        )
        v_packet.next_line_id = next_line_id

    def time_trace_vpacket_within_shell(self):
        v_packet = self.v_packet
        verysimple_numba_radial_1d_geometry = (
            self.verysimple_numba_radial_1d_geometry
        )
        verysimple_time_explosion = self.verysimple_time_explosion
        verysimple_opacity_state = self.verysimple_opacity_state
        enable_full_relativity = self.verysimple_enable_full_relativity
        continuum_processes_enabled = (
            self.verysimple_continuum_processes_enabled
        )

        # Give the vpacket a reasonable line ID
        self.v_packet_initialize_line_id(
            v_packet,
            verysimple_opacity_state,
            verysimple_time_explosion,
            enable_full_relativity,
        )

        (
            tau_trace_combined,
            distance_boundary,
            delta_shell,
        ) = vpacket.trace_vpacket_within_shell(
            v_packet,
            verysimple_numba_radial_1d_geometry,
            verysimple_time_explosion,
            verysimple_opacity_state,
            enable_full_relativity,
            continuum_processes_enabled,
        )

        assert delta_shell == 1

    def time_trace_vpacket(self):
        v_packet = self.v_packet
        verysimple_numba_radial_1d_geometry = (
            self.verysimple_numba_radial_1d_geometry
        )
        verysimple_time_explosion = self.verysimple_time_explosion
        verysimple_opacity_state = self.verysimple_opacity_state
        enable_full_relativity = self.verysimple_enable_full_relativity
        continuum_processes_enabled = (
            self.verysimple_continuum_processes_enabled
        )
        tau_russian = self.verysimple_tau_russian
        survival_probability = self.verysimple_survival_probability

        # Set seed because of RNG in trace_vpacket
        np.random.seed(1)

        # Give the vpacket a reasonable line ID
        self.v_packet_initialize_line_id(
            v_packet,
            verysimple_opacity_state,
            verysimple_time_explosion,
            enable_full_relativity,
        )

        tau_trace_combined = vpacket.trace_vpacket(
            v_packet,
            verysimple_numba_radial_1d_geometry,
            verysimple_time_explosion,
            verysimple_opacity_state,
            tau_russian,
            survival_probability,
            enable_full_relativity,
            continuum_processes_enabled,
        )

        assert v_packet.next_line_id == 2773
        assert v_packet.current_shell_id == 1

    @property
    def broken_packet(self):
        return vpacket.VPacket(
            r=1286064000000000.0,
            nu=1660428912896553.2,
            mu=0.4916053094346575,
            energy=2.474533071386993e-07,
            index=3,
            current_shell_id=0,
            next_line_id=5495,
        )

    def time_trace_bad_vpacket(self):
        broken_packet = self.broken_packet
        verysimple_numba_radial_1d_geometry = (
            self.verysimple_numba_radial_1d_geometry
        )
        enable_full_relativity = self.verysimple_enable_full_relativity
        verysimple_time_explosion = self.verysimple_time_explosion
        verysimple_opacity_state = self.verysimple_opacity_state
        continuum_processes_enabled = (
            self.verysimple_continuum_processes_enabled
        )
        tau_russian = self.verysimple_tau_russian
        survival_probability = self.verysimple_survival_probability

        vpacket.trace_vpacket(
            broken_packet,
            verysimple_numba_radial_1d_geometry,
            verysimple_time_explosion,
            verysimple_opacity_state,
            tau_russian,
            survival_probability,
            enable_full_relativity,
            continuum_processes_enabled,
        )
