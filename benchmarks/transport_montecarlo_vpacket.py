"""
Basic TARDIS Benchmark.
"""

import functools

import numpy as np

import tardis.transport.montecarlo.vpacket as vpacket
from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.frame_transformations import get_doppler_factor
from tardis.transport.montecarlo.r_packet import RPacket


class BenchmarkMontecarloMontecarloNumbaVpacket(BenchmarkBase):
    """
    Class to benchmark the single packet loop function.
    """

    repeat = 5

    def setup(self):
        self.vpacket = self.v_packet
        self.numba_radial_1d_geometry = self.verysimple_numba_radial_1d_geometry
        self.time_explosion = self.verysimple_time_explosion
        self.opacity_state = self.verysimple_opacity_state
        self.tau_russian = self.verysimple_tau_russian
        self.survival_probability = self.verysimple_survival_probability
        self.enable_full_relativity = self.verysimple_enable_full_relativity

    @functools.cached_property
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

    @functools.cached_property
    def r_packet(self):
        return RPacket(
            r=7.5e14,
            nu=4e18,
            mu=self.verysimple_packet_collection.initial_mus[0],
            energy=self.verysimple_packet_collection.initial_energies[0],
            seed=1963,
            index=0,
        )

    @functools.cache
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
        # Give the vpacket a reasonable line ID
        self.v_packet_initialize_line_id(
            self.vpacket,
            self.opacity_state,
            self.time_explosion,
            self.enable_full_relativity,
        )

        vpacket.trace_vpacket_within_shell(
            self.vpacket,
            self.numba_radial_1d_geometry,
            self.time_explosion,
            self.opacity_state,
            self.enable_full_relativity,
        )

    def time_trace_vpacket(self):
        # Set seed because of RNG in trace_vpacket
        np.random.seed(1)

        # Give the vpacket a reasonable line ID
        self.v_packet_initialize_line_id(
            self.vpacket,
            self.opacity_state,
            self.time_explosion,
            self.enable_full_relativity,
        )

        vpacket.trace_vpacket(
            self.vpacket,
            self.numba_radial_1d_geometry,
            self.time_explosion,
            self.opacity_state,
            self.tau_russian,
            self.survival_probability,
            self.enable_full_relativity,
        )

    @functools.cached_property
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

        vpacket.trace_vpacket(
            broken_packet,
            self.numba_radial_1d_geometry,
            self.time_explosion,
            self.opacity_state,
            self.tau_russian,
            self.survival_probability,
            self.enable_full_relativity,
        )

    def time_trace_vpacket_volley(self):
        vpacket.trace_vpacket_volley(
            self.r_packet,
            self.verysimple_3vpacket_collection,
            self.numba_radial_1d_geometry,
            self.time_explosion,
            self.opacity_state,
            self.enable_full_relativity,
            self.tau_russian,
            self.survival_probability,
        )
