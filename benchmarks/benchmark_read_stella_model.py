"""Basic TARDIS Benchmark."""
import os
from tardis.io.model.readers.stella import read_stella_model
from tardis.io.model import read_stella_model


class Benchmarkruntardis:
    """Class to benchmark the run_tardis function.
    """
    timeout = 200

    def setup(self):
        filename = "messa.stella.dat"
        dir_path = os.path.dirname(os.path.realpath(__file__))
        path = os.path.join(dir_path, "data", filename)
        config = read_stella_model(path)
        self.config = config

    def time_run_stella_model(self):
        stella_model = read_stella_model(self.config)
        stella_model
