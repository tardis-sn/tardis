"""
Basic TARDIS Benchmark.
"""

import numpy as np

import tardis.montecarlo.montecarlo_numba.vpacket as vpacket
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

    def v_packet_initialize_line_id(self, v_packet, opacity_state, numba_model):
        inverse_line_list_nu = opacity_state.line_list_nu[::-1]
        doppler_factor = get_doppler_factor(
            v_packet.r, v_packet.mu, numba_model.time_explosion
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
        verysimple_numba_model = self.verysimple_numba_model
        verysimple_opacity_state = self.verysimple_opacity_state

        # Give the vpacket a reasonable line ID
        self.v_packet_initialize_line_id(
            v_packet, verysimple_opacity_state, verysimple_numba_model
        )

        (
            tau_trace_combined,
            distance_boundary,
            delta_shell,
        ) = vpacket.trace_vpacket_within_shell(
            v_packet,
            verysimple_numba_radial_1d_geometry,
            verysimple_numba_model,
            verysimple_opacity_state,
        )

        assert delta_shell == 1

    def time_trace_vpacket(self):
        v_packet = self.v_packet
        verysimple_numba_radial_1d_geometry = (
            self.verysimple_numba_radial_1d_geometry
        )
        verysimple_numba_model = self.verysimple_numba_model
        verysimple_opacity_state = self.verysimple_opacity_state

        # Set seed because of RNG in trace_vpacket
        np.random.seed(1)

        # Give the vpacket a reasonable line ID
        self.v_packet_initialize_line_id(
            v_packet, verysimple_opacity_state, verysimple_numba_model
        )

        tau_trace_combined = vpacket.trace_vpacket(
            v_packet,
            verysimple_numba_radial_1d_geometry,
            verysimple_numba_model,
            verysimple_opacity_state,
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
        verysimple_numba_model = self.verysimple_numba_model
        verysimple_opacity_state = self.verysimple_opacity_state

        vpacket.trace_vpacket(
            broken_packet,
            verysimple_numba_radial_1d_geometry,
            verysimple_numba_model,
            verysimple_opacity_state,
        )
