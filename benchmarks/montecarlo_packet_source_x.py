"""
Basic TARDIS Benchmark.
"""
# import os
#
# import numpy as np
# import pandas as pd
# from asv_runner.benchmarks.mark import skip_benchmark
# from numpy.testing import assert_allclose
#
# from benchmarks.benchmark_base import BenchmarkBase
# from tardis.montecarlo import (
#     montecarlo_configuration as montecarlo_configuration,
# )
# from tardis.montecarlo.packet_source import (
#     BlackBodySimpleSource,
#     BlackBodySimpleSourceRelativistic,
# )
#
#
# # @skip_benchmark
# class BenchmarkMontecarloPacketSource(BenchmarkBase):
#     """
#     Class to benchmark the packet source function.
#     """
#
#     def __init__(self):
#         pass
#
#     @property
#     def packet_unit_test_fpath(self):
#         return f"{self.tardis_ref_path}/packet_unittest.h5"
#
#     @property
#     def blackbodysimplesource(self):
#         cls = type(self)
#         montecarlo_configuration.LEGACY_MODE_ENABLED = True
#         cls.bb = BlackBodySimpleSource(base_seed=1963, legacy_second_seed=2508)
#         yield cls.bb
#         montecarlo_configuration.LEGACY_MODE_ENABLED = False
#
#     @property
#     def blackbody_simplesource_relativistic(self, request):
#         montecarlo_configuration.LEGACY_MODE_ENABLED = True
#         bb_rel = BlackBodySimpleSourceRelativistic(
#             base_seed=1963, legacy_second_seed=2508
#         )
#         yield bb_rel
#         montecarlo_configuration.LEGACY_MODE_ENABLED = False
#
#     def time_bb_packet_sampling(self):
#         ref_df = self.tardis_ref_data["/packet_unittest/blackbody"]
#         self.bb.temperature = 10000
#         nus = self.bb.create_packet_nus(100)
#         mus = self.bb.create_packet_mus(100)
#         unif_energies = self.bb.create_packet_energies(100)
#         assert np.all(np.isclose(nus, ref_df["nus"]))
#         assert np.all(np.isclose(mus, ref_df["mus"]))
#         assert np.all(np.isclose(unif_energies, ref_df["energies"]))
#
#     # def time_bb_packet_sampling_relativistic(self):
#     #     blackbody_simplesource_relativistic.temperature = 10000
#     #     blackbody_simplesource_relativistic.beta = 0.25
#     #
#     #     nus = blackbody_simplesource_relativistic.create_packet_nus(100)
#     #     unif_energies = (
#     #         blackbody_simplesource_relativistic.create_packet_energies(100)
#     #     )
#     #     blackbody_simplesource_relativistic._reseed(2508)
#     #     mus = blackbody_simplesource_relativistic.create_packet_mus(10)
#     #
#     #     gamma = np.sqrt(1 - blackbody_simplesource_relativistic.beta ** 2) ** -1
#     #     ref_df = tardis_ref_data["/packet_unittest/blackbody"]
#     #     expected_nus = ref_df["nus"]
#     #     expected_unif_energies = ref_df["energies"] * 1.6 / gamma
#     #     expected_mus = np.array(
#     #         [
#     #             0.60420546,
#     #             0.49899691,
#     #             0.69583288,
#     #             0.96812652,
#     #             0.01544154,
#     #             0.93562304,
#     #             0.44306545,
#     #             0.77010037,
#     #             0.896973,
#     #             0.67876489,
#     #         ]
#     #     )
#     #
#     #     assert_allclose(nus, expected_nus)
#     #     assert_allclose(unif_energies, expected_unif_energies)
#     #     assert_allclose(mus, expected_mus, rtol=1e-6)
