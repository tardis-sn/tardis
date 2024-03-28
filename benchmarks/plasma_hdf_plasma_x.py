"""
Basic TARDIS Benchmark.
"""
# import numpy.testing as npt
# import pandas.testing as pdt
# from asv_runner.benchmarks.mark import skip_benchmark, parameterize
#
# from benchmarks.benchmark_base import BenchmarkBase
#
#
# @skip_benchmark
# class BenchmarkPlasmaHdfPlasma(BenchmarkBase):
#     """
#     Class to benchmark the HDF plasma function.
#     """
#
#     def __init__(self):
#         pass
#
#     ###
#     # saving and loading of plasma properties in the HDF file
#     ###
#
#     plasma_properties_list = [
#         "number_density",
#         "beta_rad",
#         "general_level_boltzmann_factor",
#         "level_boltzmann_factor",
#         "stimulated_emission_factor",
#         "t_electrons",
#         "wavelength_cm",
#         "lines_lower_level_index",
#         "ionization_data",
#         "density",
#         "atomic_mass",
#         "level_number_density",
#         "lines_upper_level_index",
#         "nu",
#         "beta_sobolev",
#         "transition_probabilities",
#         "phi",
#         "electron_densities",
#         "t_rad",
#         "selected_atoms",
#         "ion_number_density",
#         "partition_function",
#         "abundance",
#         "g_electron",
#         "g",
#         "lines",
#         "f_lu",
#         "tau_sobolevs",
#         "j_blues",
#         "metastability",
#         "w",
#         "excitation_energy",
#     ]
#
#     @parameterize({"Attributes": plasma_properties_list})
#     def time_hdf_plasma(self, attr):
#         simulation_verysimple = self.simulation_verysimple
#         # TODO: Needs to work in the class RegressionData in BenchmarkBase
#         regression_data = self.regression_data
#         if hasattr(simulation_verysimple.plasma, attr):
#             actual = getattr(simulation_verysimple.plasma, attr)
#             expected = regression_data.sync_ndarray(actual)
#             if hasattr(actual, "cgs"):
#                 actual = actual.cgs.value
#             npt.assert_allclose(actual, expected)
#
#     # def time_hdf_levels(simulation_verysimple, regression_data):
#     #     actual = simulation_verysimple.plasma.levels.to_frame()
#     #     expected = regression_data.sync_dataframe(actual)
#     #     if hasattr(actual, "cgs"):
#     #         raise ValueError("should not ever happen")
#     #     pdt.assert_frame_equal(actual, expected)
#     #
#     # SCALARS_LIST = ["time_explosion", "link_t_rad_t_electron"]
#     #
#     # @parameterize({"Attributes": SCALARS_LIST})
#     # def time_hdf_scalars(simulation_verysimple, attr, regression_data):
#     #     actual = getattr(simulation_verysimple.plasma, attr)
#     #     if hasattr(actual, "cgs"):
#     #         actual = actual.cgs.value
#     #     expected = regression_data.sync_ndarray(actual)
#     #     npt.assert_allclose(actual, expected)
#     #
#     # def time_hdf_helium_treatment(simulation_verysimple, regression_data):
#     #     actual = simulation_verysimple.plasma.helium_treatment
#     #     expected = regression_data.sync_str(actual)
#     #     assert actual == expected
#     #
#     # def time_atomic_data_uuid(simulation_verysimple, regression_data):
#     #     actual = simulation_verysimple.plasma.atomic_data.uuid1
#     #     expected = regression_data.sync_str(actual)
#     #     assert actual == expected
#     #
#     # COLLECTION_PROPERTIES = ["t_rad", "w", "density"]
#     #
#     # @parameterize({"Attributes": COLLECTION_PROPERTIES})
#     # def time_collection(simulation_verysimple, attr, regression_data):
#     #     actual = getattr(simulation_verysimple.plasma, attr)
#     #     expected = regression_data.sync_ndarray(actual)
#     #     if hasattr(actual, "cgs"):
#     #         actual = actual.cgs.value
#     #     npt.assert_allclose(actual, expected)
