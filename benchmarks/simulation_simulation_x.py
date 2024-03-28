"""
Basic TARDIS Benchmark.
"""
# from asv_runner.benchmarks.mark import parameterize
# from benchmarks.benchmark_base import BenchmarkBase
# import os
# from pathlib import Path
# import logging
# from tardis.io.configuration.config_reader import Configuration
# from tardis.simulation import Simulation
# from tardis import run_tardis
# import pandas.testing as pdt
# import numpy as np
# import pandas as pd
# import astropy.units as u
# import tardis
#
#
# # @skip_benchmark
# class BenchmarkSimulation(BenchmarkBase):
#     """
#     Class to benchmark the simulation function.
#     """
#
#     def __init__(self):
#         pass
#
#     # @property
#     def refdata(self, key):
#         print(f'key: {key}')
#         print(f'self.tardis_ref_data: {self.tardis_ref_data}')
#         return self.tardis_ref_data[f"/test_simulation/{key}"]
#         # def get_ref_data(key):
#         #     return self.tardis_ref_data[f"/test_simulation/{key}"]
#         #
#         # return get_ref_data
#
#     @property
#     def config(self):
#         return Configuration.from_yaml(
#             f"{self.example_configuration_dir}/tardis_configv1_verysimple.yml"
#         )
#
#     @property
#     def simulation_one_loop(self):
#         config = self.config
#         config.atom_data = self.atomic_data_fname
#         config.montecarlo.iterations = 2
#         config.montecarlo.no_of_packets = int(4e4)
#         config.montecarlo.last_no_of_packets = int(4e4)
#
#         simulation = Simulation.from_config(config)
#         simulation.run_convergence()
#         simulation.run_final()
#
#         print(f'1.1 simulation: {simulation}')
#         print(f'1.2 simulation.hdf_properties: {simulation.hdf_properties}')
#         print(f'1.3 simulation.simulation_state.hdf_properties: {simulation.simulation_state.hdf_properties}')
#         print(f'1.4 simulation.transport.hdf_properties: {simulation.transport.hdf_properties}')
#         # TODO: The generate reference needs to be dynamic to recognize this condition.
#         if not self.generate_reference:
#             return simulation
#         else:
#             simulation.hdf_properties = [
#                 "iterations_w",
#                 "iterations_t_rad",
#                 "iterations_electron_densities",
#                 "iterations_t_inner",
#             ]
#             simulation.simulation_state.hdf_properties = [
#                 "t_radiative",
#                 "dilution_factor",
#             ]
#             simulation.transport.hdf_properties = [
#                 "j_estimator",
#                 "nu_bar_estimator",
#                 "output_nu",
#                 "output_energy",
#             ]
#             simulation.to_hdf(
#                 self.tardis_ref_data, "", "test_simulation", overwrite=True
#             )
#             simulation.simulation_state.to_hdf(
#                 self.tardis_ref_data, "", "test_simulation", overwrite=True
#             )
#             simulation.transport.to_hdf(
#                 self.tardis_ref_data, "", "test_simulation", overwrite=True
#             )
#             raise Exception("Reference data was generated during this run.")
#
#     @parameterize({"Name": [
#             "nu_bar_estimator",
#             # "j_estimator",
#             # "t_radiative",
#             # "dilution_factor",
#             # "output_nu",
#             # "output_energy",
#         ]
#     })
#     def time_plasma_estimates(self, name):
#         print(f'self.simulation_one_loop: {self.simulation_one_loop}')
#         if name in ["nu_bar_estimator", "j_estimator"]:
#             actual = getattr(
#                 self.simulation_one_loop.transport.transport_state.radfield_mc_estimators,
#                 name,
#             )
#         elif name in ["t_radiative", "dilution_factor"]:
#             actual = getattr(self.simulation_one_loop.simulation_state, name)
#         elif name in ["output_nu", "output_energy"]:
#             OLD_TO_NEW_DICT = {
#                 "output_nu": "output_nus",
#                 "output_energy": "output_energies",
#             }
#             actual = getattr(
#                 self.simulation_one_loop.transport.transport_state.packet_collection,
#                 OLD_TO_NEW_DICT[name],
#             )
#         else:
#             try:
#                 actual = getattr(self.simulation_one_loop.transport, name)
#             except AttributeError:
#                 actual = getattr(self.simulation_one_loop.simulation_state, name)
#         if name in ["t_radiative"]:
#             # removing the quantitiness of the data
#             actual = actual.value
#         actual = pd.Series(actual)
#
#         pdt.assert_series_equal(actual, self.refdata(name), rtol=1e-5, atol=1e-8)
#
#     # @parameterize(
#     #     "name",
#     #     [
#     #         "iterations_w",
#     #         "iterations_t_rad",
#     #         "iterations_electron_densities",
#     #         "iterations_t_inner",
#     #     ],
#     # )
#     # def time_plasma_state_iterations(self, name):
#     #     actual = getattr(self.simulation_one_loop, name)
#     #
#     #     # removing the quantitiness of the data as it will screw up the comparison via pandas
#     #     if hasattr(actual, "value"):
#     #         actual = actual.value
#     #
#     #     try:
#     #         actual = pd.Series(actual)
#     #     except Exception:
#     #         actual = pd.DataFrame(actual)
#     #
#     #     if type(actual) == pd.DataFrame:
#     #         pdt.assert_frame_equal(actual, self.refdata(name), rtol=1e-5, atol=1e-8)
#     #     elif type(actual) == pd.Series:
#     #         pdt.assert_series_equal(actual, self.refdata(name))
#     #     else:
#     #         raise ValueError(f"Unknown type of actual {type(actual)}")
#     #
#     # @pytest.fixture(scope="module")
#     # def simulation_without_loop(atomic_data_fname, config):
#     #     config.atom_data = atomic_data_fname
#     #     config.montecarlo.iterations = 2
#     #     return Simulation.from_config(config)
#     #
#     # def time_plasma_state_storer_store(
#     #         atomic_data_fname, config, simulation_without_loop
#     # ):
#     #     simulation = simulation_without_loop
#     #
#     #     w_test = np.linspace(0, 1, 20)
#     #     t_rad_test = np.linspace(12000, 9000, 20) * u.K
#     #     electron_densities_test = pd.Series(np.linspace(1e7, 1e6, 20))
#     #     t_inner_test = 12500 * u.K
#     #
#     #     simulation.store_plasma_state(
#     #         1, w_test, t_rad_test, electron_densities_test, t_inner_test
#     #     )
#     #
#     #     np.testing.assert_allclose(simulation.iterations_w[1, :], w_test)
#     #     np.testing.assert_allclose(simulation.iterations_t_rad[1, :], t_rad_test)
#     #     np.testing.assert_allclose(
#     #         simulation.iterations_electron_densities[1, :], electron_densities_test
#     #     )
#     #     np.testing.assert_allclose(simulation.iterations_t_inner[1], t_inner_test)
#     #
#     # def time_plasma_state_storer_reshape(
#     #         atomic_data_fname, config, simulation_without_loop
#     # ):
#     #     simulation = simulation_without_loop
#     #     simulation.reshape_plasma_state_store(0)
#     #
#     #     assert simulation.iterations_t_rad.shape == (1, 20)
#     #     assert simulation.iterations_w.shape == (1, 20)
#     #     assert simulation.iterations_electron_densities.shape == (1, 20)
#     #     assert simulation.iterations_t_inner.shape == (1,)
#     #
#     # def time_version_tag(simulation_without_loop):
#     #     simulation = simulation_without_loop
#     #     assert simulation.version == tardis.__version__
