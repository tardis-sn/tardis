"""
Basic TARDIS Benchmark.
"""
from astropy import units as u
from asv_runner.benchmarks.mark import skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.io.configuration.config_reader import ConfigurationNameSpace


# @skip_benchmark
class BenchmarkIoConfigurationConfigurationNamespace(BenchmarkBase):
    """
    Class to benchmark the run_tardis function.
    """

    def __init__(self):
        self.simple_config_dict = {
            "a": {"b": {"param1": 1, "param2": [0, 1, 2 * u.km], "param3": 4.0 * u.km}}
        }

    @property
    def config_ns(self):
        """Return example instance of `ConfigurationNameSpace` class."""
        config_name_space = ConfigurationNameSpace(self.simple_config_dict)
        return config_name_space

    def time_simple_configuration_namespace(self):
        self.config_ns.a.b.param1 = 2
        self.config_ns["a"]["b"]["param1"] = 3

    def time_quantity_configuration_namespace(self):
        self.config_ns.a.b.param3 = 3
        self.config_ns.a.b.param3 = 5000 * u.m

    def time_access_with_config_item_string(self):
        self.config_ns.set_config_item("a.b.param1", 2)

    def time_set_with_config_item_string_quantity(self):
        self.config_ns.set_config_item("a.b.param3", 2)

    def time_get_with_config_item_string_item_access(self):
        self.config_ns.get_config_item("a.b.param2.item0")
        self.config_ns.get_config_item("a.b.param2.item1")

    def time_set_with_config_item_string_item_access(self):
        self.config_ns.set_config_item("a.b.param2.item0", 2)
        self.config_ns.get_config_item("a.b.param2.item0")

    def time_set_with_config_item_string_item_access_quantity(self):
        self.config_ns.set_config_item("a.b.param2.item2", 7)
        self.config_ns.get_config_item("a.b.param2.item2")

    def time_config_namespace_copy(self):
        config_ns2 = self.config_ns.deepcopy()
        config_ns2.a.b.param1 = 2
