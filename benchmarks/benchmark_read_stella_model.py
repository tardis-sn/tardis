"""Basic TARDIS Benchmark."""
import os
from tardis.io.model import read_stella_model


class Benchmark_read_stella_model:
    """Class to benchmark the read_stella_model function.
    """
    timeout = 200

    def setup(self):
        filename = "messa.stella.dat"
        dir_path = os.path.dirname(os.path.realpath(__file__))
        path = os.path.join(dir_path, "data", filename)
        self.path = path

    def time_run_stella_model(self):
        stella_model = read_stella_model(self.path)
        stella_model
