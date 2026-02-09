import functools
from os import environ as env
from pathlib import Path

from tardis.io.configuration.config_reader import Configuration


class BenchmarkBase:
    # the total time for all the repetitions for each benchmark.
    timeout = 600

    @staticmethod
    def get_relative_path(partial_path: str) -> str:
        path = Path(__file__).resolve().parent

        for target in Path(partial_path).parts:
            path = path / target

        return str(path)

    def get_absolute_path(self, partial_path):
        partial_path = "../" + partial_path

        return self.get_relative_path(partial_path)

    @functools.cached_property
    def tardis_ref_path(self):
        ref_data_path = Path(
            Path(__file__).parent.parent, env.get("TARDIS_REF_PATH")
        ).resolve()
        return ref_data_path

    @functools.cached_property
    def atomic_data_fname(self):
        atomic_data_fname = (
            f"{self.tardis_ref_path}/kurucz_cd23_chianti_H_He_latest.h5"
        )

        if not Path(atomic_data_fname).exists():
            atom_data_missing_str = (
                f"{atomic_data_fname} atomic datafiles does not seem to exist"
            )
            raise Exception(atom_data_missing_str)

        return atomic_data_fname

    @functools.cached_property
    def example_configuration_dir(self):
        return self.get_absolute_path("tardis/io/configuration/tests/data")

    @functools.cached_property
    def config_verysimple(self):
        return Configuration.from_yaml(
            f"{self.example_configuration_dir}/tardis_configv1_verysimple.yml"
        )

    @functools.cached_property
    def config_rpacket_tracking(self):
        config = Configuration.from_yaml(
            f"{self.example_configuration_dir}/tardis_configv1_verysimple.yml"
        )
        config.montecarlo.tracking.track_rpacket = True
        return config
