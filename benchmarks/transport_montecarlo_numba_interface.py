"""
Basic TARDIS Benchmark.
"""

import numpy as np
from asv_runner.benchmarks.mark import parameterize

import tardis.opacities.opacity_state as numba_interface
from benchmarks.benchmark_base import BenchmarkBase


class BenchmarkMontecarloMontecarloNumbaNumbaInterface(BenchmarkBase):
    """
    Class to benchmark the numba interface function.
    """

    @parameterize({"Input params": ["scatter", "macroatom", "downbranch"]})
    def time_opacity_state_initialize(self, input_params):
        line_interaction_type = input_params
        plasma = self.nb_simulation_verysimple.plasma
        numba_interface.opacity_state_initialize(
            plasma,
            line_interaction_type,
            self.verysimple_disable_line_scattering,
            self.verysimple_continuum_processes_enabled,
        )

        if line_interaction_type == "scatter":
            np.zeros(1, dtype=np.int64)

    def time_VPacketCollection_add_packet(self):
        verysimple_3vpacket_collection = self.verysimple_3vpacket_collection
        assert verysimple_3vpacket_collection.length == 0

        nus = [3.0e15, 0.0, 1e15, 1e5]
        energies = [0.4, 0.1, 0.6, 1e10]
        initial_mus = [0.1, 0, 1, 0.9]
        initial_rs = [3e42, 4.5e45, 0, 9.0e40]
        last_interaction_in_nus = np.array(
            [3.0e15, 0.0, 1e15, 1e5], dtype=np.float64
        )
        last_interaction_types = np.array([1, 1, 3, 2], dtype=np.int64)
        last_interaction_in_ids = np.array([100, 0, 1, 1000], dtype=np.int64)
        last_interaction_out_ids = np.array(
            [1201, 123, 545, 1232], dtype=np.int64
        )
        last_interaction_shell_ids = np.array([2, -1, 6, 0], dtype=np.int64)

        for (
            nu,
            energy,
            initial_mu,
            initial_r,
            last_interaction_in_nu,
            last_interaction_type,
            last_interaction_in_id,
            last_interaction_out_id,
            last_interaction_shell_id,
        ) in zip(
            nus,
            energies,
            initial_mus,
            initial_rs,
            last_interaction_in_nus,
            last_interaction_types,
            last_interaction_in_ids,
            last_interaction_out_ids,
            last_interaction_shell_ids,
        ):
            verysimple_3vpacket_collection.add_packet(
                nu,
                energy,
                initial_mu,
                initial_r,
                last_interaction_in_nu,
                last_interaction_type,
                last_interaction_in_id,
                last_interaction_out_id,
                last_interaction_shell_id,
            )

        assert verysimple_3vpacket_collection.length == 9
