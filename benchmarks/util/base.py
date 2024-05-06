from os.path import dirname, realpath, join
from pathlib import Path, PosixPath


class Base:
    @staticmethod
    def get_path(partial_path: str) -> Path:
        base_path = dirname(realpath(__file__))
        path = Path(base_path) / Path(partial_path)
        return path

    @property
    def tardis_ref_path(self) -> Path:
        # TODO: This route is fixed but needs to get from the arguments given in the command line.
        #       /app/tardis-refdata
        return Path("/app/tardis-refdata")

    @property
    def example_configuration_dir(self) -> Path:
        return self.get_path("../../tardis/io/configuration/tests/data")
