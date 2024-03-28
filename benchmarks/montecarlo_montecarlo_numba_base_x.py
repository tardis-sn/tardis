"""
Basic TARDIS Benchmark.
"""
# from copy import deepcopy
#
# from asv_runner.benchmarks.mark import skip_benchmark
#
# from benchmarks.benchmark_base import BenchmarkBase
# from tardis.montecarlo import (
#     montecarlo_configuration as montecarlo_configuration,
# )
# from tardis.simulation import Simulation
#
#
# # @skip_benchmark
# class BenchmarkMontecarloMontecarloNumbaBase(BenchmarkBase):
#     """
#     Class to benchmark the montecarlo numba base function.
#     """
#
#     def __init__(self):
#         pass
#
#     @property
#     def montecarlo_main_loop_config(self):
#         montecarlo_configuration.LEGACY_MODE_ENABLED = True
#         # Setup model config from verysimple
#
#         config_montecarlo_1e5_verysimple = self.config_verysimple
#         config_montecarlo_1e5_verysimple.montecarlo.last_no_of_packets = 1e5
#         config_montecarlo_1e5_verysimple.montecarlo.no_of_virtual_packets = 0
#         config_montecarlo_1e5_verysimple.montecarlo.iterations = 1
#         config_montecarlo_1e5_verysimple.plasma.line_interaction_type = "macroatom"
#
#         del config_montecarlo_1e5_verysimple["config_dirname"]
#         return config_montecarlo_1e5_verysimple
#
#     def time_montecarlo_main_loop(self):
#         atomic_dataset = deepcopy(self.atomic_dataset)
#         montecarlo_main_loop_simulation = Simulation.from_config(
#             self.montecarlo_main_loop_config,
#             atom_data=atomic_dataset,
#             virtual_packet_logging=False,
#         )
#         montecarlo_main_loop_simulation.run_convergence()
#         montecarlo_main_loop_simulation.run_final()
#
#         expected_hdf_store = self.regression_data.sync_hdf_store(
#             montecarlo_main_loop_simulation
#         )
#
#         expected_hdf_store.close()
#
#     def time_montecarlo_main_loop_vpacket_log(self):
#         atomic_dataset = deepcopy(self.atomic_dataset)
#         self.montecarlo_main_loop_config.montecarlo.no_of_virtual_packets = 5
#
#         montecarlo_main_loop_simulation = Simulation.from_config(
#             self.montecarlo_main_loop_config,
#             atom_data=atomic_dataset,
#             virtual_packet_logging=True,
#         )
#         montecarlo_main_loop_simulation.run_convergence()
#         montecarlo_main_loop_simulation.run_final()
#
#         assert montecarlo_configuration.ENABLE_VPACKET_TRACKING == True
#
#         expected_hdf_store = self.regression_data.sync_hdf_store(
#             montecarlo_main_loop_simulation
#         )
#
#         expected_hdf_store.close()
