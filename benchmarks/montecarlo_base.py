"""
Basic TARDIS Benchmark.
"""
import pandas as pd
from asv_runner.benchmarks.mark import skip_benchmark, parameterize
from numpy.testing import assert_almost_equal

from benchmarks.benchmark_base import BenchmarkBase


# @skip_benchmark
class BenchmarkMontecarloBase(BenchmarkBase):
    """
    Class to benchmark the montecarlo base function.
    """

    def __init__(self):
        pass

    transport_state_properties = [
        "output_nu",
        "output_energy",
        "nu_bar_estimator",
        "j_estimator",
        "montecarlo_virtual_luminosity",
        "packet_luminosity",
        # These are nested properties that should be tested differently
        # "spectrum",
        # "spectrum_virtual",
        # "spectrum_reabsorbed",
        # This is a scalar and should be tested differently
        # "time_of_simulation",
        "emitted_packet_mask",
        "last_interaction_type",
        "last_interaction_in_nu",
        "last_line_interaction_out_id",
        "last_line_interaction_in_id",
        "last_line_interaction_shell_id",
        "virt_packet_nus",
        "virt_packet_energies",
        "virt_packet_initial_rs",
        "virt_packet_initial_mus",
        "virt_packet_last_interaction_in_nu",
        "virt_packet_last_interaction_type",
        "virt_packet_last_line_interaction_in_id",
        "virt_packet_last_line_interaction_out_id",
        "virt_packet_last_line_interaction_shell_id",
    ]

    @parameterize({"Transport state properties": transport_state_properties})
    def time_hdf_transport_state(self, attr):
        hdf_file_path = self.hdf_file_path
        self.simulation_verysimple_vpacket_tracking.transport.to_hdf(
            hdf_file_path, name="transport", overwrite=True
        )
        self.simulation_verysimple_vpacket_tracking.transport.transport_state.to_hdf(
            hdf_file_path, name="transport_state", overwrite=True
        )
        actual = getattr(
            self.simulation_verysimple_vpacket_tracking.transport.transport_state, attr
        )
        if hasattr(actual, "cgs"):
            actual = actual.cgs.value
        path = f"transport_state/{attr}"
        expected = pd.read_hdf(hdf_file_path, path)
        assert_almost_equal(actual, expected.values)
