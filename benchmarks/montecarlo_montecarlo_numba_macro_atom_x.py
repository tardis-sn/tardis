"""
Basic TARDIS Benchmark.
"""
# from asv_runner.benchmarks.mark import parameterize, skip_benchmark
#
# import tardis.montecarlo.montecarlo_numba.macro_atom as macro_atom
# from benchmarks.benchmark_base import BenchmarkBase
#
#
# @skip_benchmark
# class BenchmarkMontecarloMontecarloNumbaMacroAtom(BenchmarkBase):
#     """
#     Class to benchmark the macro atom function.
#     """
#
#     def __init__(self):
#         pass
#
#     @parameterize({
#         "Seed": [
#             1963, 1, 2111963, 10000
#         ],
#         "Expected": [
#             10015, 9993, 17296, 9993
#         ]
#     })
#     def time_macro_atom(self, seed, expected):
#         # TODO: Check how to set the seed outside the PyTest.
#         self.set_seed_fixture(seed)
#         self.static_packet.initialize_line_id(
#             self.verysimple_opacity_state, self.verysimple_numba_model
#         )
#         activation_level_id = self.verysimple_opacity_state.line2macro_level_upper[
#             self.static_packet.next_line_id
#         ]
#         result, transition_type = macro_atom.macro_atom(
#             activation_level_id,
#             self.static_packet.current_shell_id,
#             self.verysimple_opacity_state,
#         )
#         assert result == expected
#         assert transition_type == -1  # line transition
