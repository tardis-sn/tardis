"""
Basic TARDIS Benchmark.
"""
import json
import os

import numpy as np
from asv_runner.benchmarks.mark import skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.io.model.readers import arepo


# @skip_benchmark
class BenchmarkIoModelReadersArepoParser(BenchmarkBase):
    """
    Class to benchmark the arepo parser function.
    """

    def __init__(self):
        self.test_filename = self.create_temporal_file('arepo_parser_test.csvy')

    @property
    def arepo_snapshot_fname(self):
        return f"{self.tardis_ref_path}/arepo_data/arepo_snapshot.json"

    def get_cone_csvy_model(self):
        with open(self.arepo_snapshot_fname, "r") as json_file:
            data = json.loads(json.load(json_file))

        pos, vel, rho, mass, nucs, time = (
            data["pos"],
            data["vel"],
            data["rho"],
            data["mass"],
            data["xnuc"],
            data["time"],
        )
        pos = np.array(pos)
        vel = np.array(vel)
        rho = np.array(rho)
        mass = np.array(mass)

        # The nuclear data should be in a dict where each element has its own entry (with the key being the element name)
        xnuc = {
            "ni56": np.array(nucs[0]),
            "si28": np.array(nucs[1]),
        }

        profile = arepo.ConeProfile(pos, vel, rho, mass, xnuc, time)

        profile.create_profile(
            opening_angle=40, inner_radius=1e11, outer_radius=2e11
        )

        testfile = profile.export(
            20, f"{self.test_filename}"
        )

        with open(testfile, "r") as file:
            data = file.readlines()

        print(testfile)
        os.remove(testfile)

        return data

    def get_full_csvy_model(self):
        with open(self.arepo_snapshot_fname, "r") as json_file:
            data = json.loads(json.load(json_file))

        pos, vel, rho, mass, nucs, time = (
            data["pos"],
            data["vel"],
            data["rho"],
            data["mass"],
            data["xnuc"],
            data["time"],
        )
        pos = np.array(pos)
        vel = np.array(vel)
        rho = np.array(rho)
        mass = np.array(mass)

        # The nuclear data should be in a dict where each element has its own entry (with the key being the element name)
        xnuc = {
            "ni56": np.array(nucs[0]),
            "si28": np.array(nucs[1]),
        }

        profile = arepo.FullProfile(pos, vel, rho, mass, xnuc, time)

        profile.create_profile(inner_radius=1e11, outer_radius=2e11)

        testfile = profile.export(
            20, f"{self.test_filename}"
        )

        with open(testfile, "r") as file:
            data = file.readlines()

        os.remove(testfile)

        return data

    def time_cone_profile(self):
        self.get_cone_csvy_model()

    def time_full_profile(self):
        self.get_full_csvy_model()
