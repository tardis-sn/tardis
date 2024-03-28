"""
Basic TARDIS Benchmark.
"""
from asv_runner.benchmarks.mark import skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.io.configuration.config_validator import validate_dict
from tardis.io.model.readers import csvy


# @skip_benchmark
class BenchmarkIoModelReadersCsvyReader(BenchmarkBase):
    """
    Class to benchmark the csvy reader function.
    """

    def __init__(self):
        pass

    @property
    def csvy_full_fname(self):
        return f"{self.example_model_file_dir}/csvy_full.csvy"

    @property
    def csvy_nocsv_fname(self):
        return f"{self.example_model_file_dir}/csvy_nocsv.csvy"

    @property
    def csvy_missing_fname(self):
        return f"{self.example_model_file_dir}/csvy_missing.csvy"

    def time_csvy_finds_csv_first_line(self):
        csvy.load_csvy(self.csvy_full_fname)

    def time_csv_colnames_equiv_datatype_fields(self):
        yaml_dict, csv = csvy.load_csvy(self.csvy_full_fname)
        datatype_names = [od["name"] for od in yaml_dict["datatype"]["fields"]]
        for key in csv.columns:
            assert key in datatype_names
        for name in datatype_names:
            assert name in csv.columns

    def time_csvy_nocsv_data_is_none(self):
        csvy.load_csvy(self.csvy_nocsv_fname)

    def time_csvy_full_with_validate_dictionary(self):
        yaml_dict, _ = csvy.load_csvy(self.csvy_full_fname)
        validate_dict(
            yaml_dict,
            schemapath=self.get_absolute_path('tardis/io/configuration/schemas/csvy_model.yml'),
        )
