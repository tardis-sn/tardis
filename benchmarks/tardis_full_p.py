"""
Basic TARDIS Benchmark.
"""
import numpy.testing as npt
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from asv_runner.benchmarks.mark import skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis import run_tardis
from tardis.io.configuration.config_reader import Configuration
from tardis.simulation.base import Simulation


# @skip_benchmark
class BenchmarkTardisFull(BenchmarkBase):
    """
    Class to benchmark the TARDIS full function.
    """

    def __init__(self):
        pass

    # @skip_benchmark
    def time_run_tardis_from_config_obj(self):
        """
        Tests whether the run_tardis function can take in the Configuration object
        as arguments
        """
        config = Configuration.from_yaml(
            f"{self.example_configuration_dir}/tardis_configv1_verysimple.yml"
        )
        config["atom_data"] = self.atomic_data_fname

        try:
            run_tardis(config)
        except Exception as e:
            raise Exception(str(e.args[0]))


# # @skip_benchmark
# class BenchmarkTardisFullTransportSimple(BenchmarkBase):
#     """
#     Very simple run
#     """
#
#     def __init__(self):
#         self.name = "test_transport_simple"
#
#     @property
#     def transport(self):
#         config = Configuration.from_yaml(
#             f"{self.example_configuration_dir}/tardis_configv1_verysimple.yml"
#         )
#         config["atom_data"] = self.atomic_data_fname
#
#         simulation = Simulation.from_config(config)
#         simulation.run_convergence()
#         simulation.run_final()
#
#         if not self.generate_reference:
#             return simulation.transport
#         else:
#             simulation.transport.hdf_properties = [
#                 "j_blue_estimator",
#                 "spectrum",
#                 "spectrum_virtual",
#             ]
#             simulation.transport.to_hdf(
#                 self.tardis_ref_data, "", self.name, overwrite=True
#             )
#             raise Exception("Reference data was generated during this run.")
#
#     def refdata(self, key):
#         # TODO: Fails because the ref data problem.
#         #        File "/app/code/benchmarks/tardis_full.py", line 75, in refdata
#         #                    return self.tardis_ref_data[f"{self.name}/{key}"]
#         #                           ~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^
#         #        TypeError: 'generator' object is not subscriptable
#         return self.tardis_ref_data[f"{self.name}/{key}"]
#
#     def time_j_blue_estimators(self):
#         j_blue_estimator = self.refdata("j_blue_estimator").values
#
#         npt.assert_allclose(
#             self.transport.transport_state.radfield_mc_estimators.j_blue_estimator,
#             j_blue_estimator,
#         )
#
#     def time_spectrum(self):
#         luminosity = u.Quantity(self.refdata("spectrum/luminosity"), "erg /s")
#
#         assert_quantity_allclose(
#             self.transport.transport_state.spectrum.luminosity, luminosity
#         )
#
#     def time_virtual_spectrum(self):
#         luminosity = u.Quantity(
#             self.refdata("spectrum_virtual/luminosity"), "erg /s"
#         )
#
#         assert_quantity_allclose(
#             self.transport.transport_state.spectrum_virtual.luminosity, luminosity
#         )
