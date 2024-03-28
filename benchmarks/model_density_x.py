"""
Basic TARDIS Benchmark.
"""
# import pandas as pd
# from asv_runner.benchmarks.mark import skip_benchmark
#
# from benchmarks.benchmark_base import BenchmarkBase
#
#
# # @skip_benchmark
# class BenchmarkModelDensity(BenchmarkBase):
#     """
#     Class to benchmark the density function.
#     """
#
#     def __init__(self):
#         pass
#
#     # @property
#     # def to_hdf_buffer(self):
#     #     self.simulation_verysimple.simulation_state.to_hdf(self.hdf_file_path, overwrite=True)
#
#     def time_hdf_density_0(self):
#         path = "simulation_state/density"
#         pd.read_hdf(self.hdf_file_path, path)
#
#     def time_hdf_time_0(self):
#         path = "simulation_state/scalars"
#         pd.read_hdf(self.hdf_file_path, path)
