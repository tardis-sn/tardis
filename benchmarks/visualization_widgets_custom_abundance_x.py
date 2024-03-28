"""
Basic TARDIS Benchmark.
"""
# from asv_runner.benchmarks.mark import skip_benchmark, parameterize
# from benchmarks.benchmark_base import BenchmarkBase
# from pathlib import Path
# import tardis
# import numpy as np
# import numpy.testing as npt
# from tardis.tests.test_util import monkeysession
# from tardis.visualization.widgets.custom_abundance import (
#     CustomAbundanceWidgetData,
#     CustomYAML,
#     CustomAbundanceWidget,
# )
#
#
# @skip_benchmark
# class BenchmarkVisualizationWidgetsCustomAbundance(BenchmarkBase):
#     """
#     Class to benchmark the custom abundance function.
#     """
#
#     def __init__(self):
#         pass
#
#     @property
#     def yml_data(self):
#         yml_path = f"{self.example_configuration_dir}/tardis_configv1_verysimple.yml"
#
#         return CustomAbundanceWidgetData.from_yml(
#             yml_path, atom_data=self.atomic_dataset
#         )
#
#     # @property
#     # def caw(yml_data, monkeysession):
#     #     # TODO: Check how to generate monkey session in the benchmarks.
#     #     caw = CustomAbundanceWidget(yml_data)
#     #     monkeysession.setattr(
#     #         "tardis.visualization.widgets.custom_abundance.is_notebook",
#     #         lambda: True,
#     #     )
#     #     caw.display()
#     #     return caw
#
#
# # class BenchmarkCustomAbundanceWidgetData(BenchmarkVisualizationWidgetsCustomAbundance):
# #     def time_get_symbols(self):
# #         symbols = self.yml_data.get_symbols()
# #         npt.assert_array_equal(symbols, ["O", "Mg", "Si", "S", "Ar", "Ca"])
#
#
# # class BenchmarkCustomAbundanceWidget(BenchmarkVisualizationWidgetsCustomAbundance):
# #     def time_update_input_item_value(self):
# #         caw = self.caw
# #         caw.update_input_item_value(0, 0.33333)
# #         assert caw.input_items[0].value == 0.333
# #
# #     def time_read_abundance(self):
# #         caw = self.caw
# #         caw.data.abundance[0] = 0.2
# #         caw.read_abundance()
# #         for i in range(caw.no_of_elements):
# #             assert caw.input_items[i].value == 0.2
# #
# #     def time_update_abundance_plot(self):
# #         caw = self.caw
# #         caw.data.abundance.iloc[0, :] = 0.2
# #         caw.update_abundance_plot(0)
# #
# #         npt.assert_array_equal(
# #             caw.fig.data[2].y, np.array([0.2] * (caw.no_of_shells + 1))
# #         )
# #
# #     def time_bound_locked_sum_to_1(self):
# #         caw = self.caw
# #         # bound checked input to 1
# #         caw.checks[0].value = True
# #         caw.checks[1].value = True
# #         caw.input_items[0].value = 0.5
# #         caw.input_items[1].value = 0.6
# #         assert caw.input_items[1].value == 0.5
# #
# #         # bound to 1 when input is checked
# #         caw.checks[2].value = True
# #         assert caw.input_items[2].value == 0
# #
# #     @parameterize({"Parameters": [
# #         {
# #             "v0": 11000,
# #             "v1": 11450,
# #             "expected": "hidden"
# #         },
# #         {
# #             "v0": 11100,
# #             "v1": 11200,
# #             "expected": "hidden"
# #         },
# #         {
# #             "v0": 11000,
# #             "v1": 11451,
# #             "expected": "visible"
# #         },
# #     ]})
# #     def time_overwrite_existing_shells(self, parameters):
# #         caw = self.caw
# #         v0 = parameters['v0']
# #         v1 = parameters['v1']
# #         expected = parameters['expected']
# #         caw.input_v_start.value = v0
# #         caw.input_v_end.value = v1
# #
# #         assert caw.overwrite_warning.layout.visibility == expected
# #
# #     @parameterize({"Parameters": [
# #         {
# #             "multishell_edit": False,
# #             "inputs": [0, 0, 0, 0, 0, 0, 0],
# #             "locks": [False] * 6,
# #             "expected": [0, 0, 0, 0, 0, 0, 0]
# #         },
# #         {
# #             "multishell_edit": False,
# #             "inputs": [0.1, 0.2, 0, 0, 0, 0, 0],
# #             "locks": [True] + [False] * 5,
# #             "expected": [0.1, 0.9, 0, 0, 0, 0, 0]
# #         },
# #         {
# #             "multishell_edit": False,
# #             "inputs": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
# #             "locks": [False] * 6,
# #             "expected": [0.0476, 0.0952, 0.143, 0.19, 0.238, 0.286]
# #         },
# #         {
# #             "multishell_edit": False,
# #             "inputs": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
# #             "locks": [True] * 2 + [False] * 4,
# #             "expected": [0.1, 0.2, 0.117, 0.156, 0.194, 0.233]
# #         },
# #         {
# #             "multishell_edit": True,
# #             "inputs": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
# #             "locks": [True] * 2 + [False] * 4,
# #             "expected": [0.1, 0.2, 0.117, 0.156, 0.194, 0.233]
# #         },
# #     ]})
# #     def time_on_btn_norm(self, parameters):
# #         caw = self.caw
# #         multishell_edit = parameters['multishell_edit']
# #         inputs = parameters['inputs']
# #         locks = parameters['locks']
# #         expected = parameters['expected']
# #
# #         if multishell_edit:
# #             caw.rbs_multi_apply.index = 0
# #             for i, item in enumerate(caw.input_items):
# #                 item.value = inputs[i]
# #                 caw.checks[i].value = locks[i]
# #
# #             caw.on_btn_norm(None)
# #
# #             for i, item in enumerate(caw.input_items):
# #                 assert item.value == expected[i]
# #
# #             start_no = caw.irs_shell_range.value[0]
# #             end_no = caw.irs_shell_range.value[1]
# #
# #             for i, v in enumerate(expected):
# #                 line = caw.fig.data[2 + i].y[start_no - 1: end_no]
# #                 unique_v = set(line)
# #                 assert len(unique_v) == 1
# #                 unique_v = float("{:.3g}".format(list(unique_v)[0]))
# #                 assert unique_v == v
# #         else:
# #             for i, item in enumerate(caw.input_items):
# #                 item.value = inputs[i]
# #                 caw.checks[i].value = locks[i]
# #
# #             caw.on_btn_norm(None)
# #
# #             for i, item in enumerate(caw.input_items):
# #                 assert item.value == expected[i]
#
#
# # class BenchmarkCustomYAML(BenchmarkVisualizationWidgetsCustomAbundance):
# #     def time_create_fields_dict(self):
# #         custom_yaml = CustomYAML("test", 0, 0, 0, 0)
# #         custom_yaml.create_fields_dict(["H", "He"])
# #         datatype_dict = {
# #             "fields": [
# #                 {"name": "velocity", "unit": "km/s"},
# #                 {"name": "density", "unit": "g/cm^3"},
# #                 {"name": "H", "desc": "fractional H abundance"},
# #                 {"name": "He", "desc": "fractional He abundance"},
# #             ]
# #         }
# #
# #         assert custom_yaml.datatype == datatype_dict
