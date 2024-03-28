"""
Basic TARDIS Benchmark.
"""
import astropy.units as u
from asv_runner.benchmarks.mark import parameterize, skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.io.configuration.config_reader import Configuration
from tardis.simulation import Simulation


# @skip_benchmark
class BenchmarkPlasmaPlasmaVBoundary(BenchmarkBase):
    """
    Class to benchmark the plasma vBoundary function.
    """

    def __init__(self):
        self.DATA_PATH = self.get_absolute_path("tardis/plasma/tests/data")

    @property
    def config_init_trad_fname(self):
        return f"{self.DATA_PATH}/config_init_trad.yml"

    @parameterize({"Parameters": [
        {
            "v_inner_boundary": 3350,
            "v_outer_boundary": 3650
         },
        {
            "v_inner_boundary": 2900,
            "v_outer_boundary": 3750
         },
        {
            "v_inner_boundary": 2900,
            "v_outer_boundary": 3850
         },
        {
            "v_inner_boundary": 2900,
            "v_outer_boundary": 3900
         },
        {
            "v_inner_boundary": 2950,
            "v_outer_boundary": 3750
         },
        {
            "v_inner_boundary": 2950,
            "v_outer_boundary": 3850
         },
        {
            "v_inner_boundary": 2950,
            "v_outer_boundary": 3900
         },
        {
            "v_inner_boundary": 3050,
            "v_outer_boundary": 3750
         },
        {
            "v_inner_boundary": 3050,
            "v_outer_boundary": 3850
         },
        {
            "v_inner_boundary": 3050,
            "v_outer_boundary": 3900
         },
        {
            "v_inner_boundary": 3150,
            "v_outer_boundary": 3750
         },
        {
            "v_inner_boundary": 3150,
            "v_outer_boundary": 3850
         },
        {
            "v_inner_boundary": 3150,
            "v_outer_boundary": 3900
         },
    ]})
    def time_plasma_vboundary(self, parameters):
        v_inner_boundary = parameters['v_inner_boundary']
        v_outer_boundary = parameters['v_outer_boundary']
        tardis_config = Configuration.from_yaml(self.config_init_trad_fname)
        tardis_config.atom_data = self.atomic_data_fname
        tardis_config.model.structure.v_inner_boundary = (
                v_inner_boundary * u.km / u.s
        )
        tardis_config.model.structure.v_outer_boundary = (
                v_outer_boundary * u.km / u.s
        )
        Simulation.from_config(tardis_config)
