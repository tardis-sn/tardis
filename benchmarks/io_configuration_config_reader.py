"""
Basic TARDIS Benchmark.
"""
import pandas as pd
from astropy.units import Quantity
from asv_runner.benchmarks.mark import parameterize, skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.io.configuration import config_reader
from tardis.io.configuration.config_reader import Configuration
from tardis.plasma.standard_plasmas import assemble_plasma


# @skip_benchmark
class BenchmarkIoConfigurationConfigReader(BenchmarkBase):
    """
    Class to benchmark the config reader function.
    """

    @staticmethod
    def time_convergence_section_parser():
        convergence_section = {
            "type": "damped",
            "lock_t_inner_cyles": 1,
            "t_inner_update_exponent": -0.5,
            "damping_constant": 0.5,
            "threshold": 0.05,
            "fraction": 0.8,
            "hold_iterations": 3,
            "t_rad": {"damping_constant": 1.0},
        }

        config_reader.Configuration.parse_convergence_section(
            convergence_section
        )

    def time_from_config_dict(self):
        Configuration.from_config_dict(
            self.tardis_config_verysimple, validate=True, config_dirname="test"
        )

    def time_config_hdf(self):
        # TODO: Delete this temporal file after ASV runs the benchmarks.
        hdf_file_path = self.hdf_file_path
        expected = Configuration.from_config_dict(
            self.tardis_config_verysimple, validate=True, config_dirname="test"
        )
        expected.to_hdf(hdf_file_path, overwrite=True)
        pd.read_hdf(hdf_file_path, key="/simulation/config")

    def time_model_section_config(self):
        Configuration.from_config_dict(
            self.tardis_config_verysimple, validate=True, config_dirname="test"
        )

        self.tardis_config_verysimple["model"]["structure"]["velocity"][
            "start"
        ] = Quantity("2.0e4 km/s")
        self.tardis_config_verysimple["model"]["structure"]["velocity"][
            "stop"
        ] = Quantity("1.1e4 km/s")

        Configuration.from_config_dict(
            self.tardis_config_verysimple, validate=True, config_dirname="test"
        )

    def time_supernova_section_config(self):
        self.tardis_config_verysimple["supernova"]["time_explosion"] = Quantity(
            "-10 day"
        )
        Configuration.from_config_dict(
            self.tardis_config_verysimple, validate=True, config_dirname="test"
        )

        self.tardis_config_verysimple["supernova"]["time_explosion"] = Quantity("10 day")
        self.tardis_config_verysimple["supernova"][
            "luminosity_wavelength_start"
        ] = Quantity("15 angstrom")
        self.tardis_config_verysimple["supernova"][
            "luminosity_wavelength_end"
        ] = Quantity("0 angstrom")
        Configuration.from_config_dict(
            self.tardis_config_verysimple, validate=True, config_dirname="test"
        )

    @parameterize({"Plasma values": ("initial_t_inner", "initial_t_rad")})
    def time_plasma_section_config(self, key):
        self.tardis_config_verysimple["plasma"][key] = Quantity("-100 K")
        Configuration.from_config_dict(
            self.tardis_config_verysimple, validate=True, config_dirname="test"
        )

    def time_plasma_nlte_section_root_config(self):
        self.nlte.tardis_config_verysimple_nlte["plasma"]["continuum_interaction"]["species"] = ["He I", ]
        self.nlte.tardis_config_verysimple_nlte["plasma"]["nlte_ionization_species"] = ["H I"]
        config = Configuration.from_config_dict(self.nlte.tardis_config_verysimple_nlte)
        assemble_plasma(config, self.nlte.nlte_raw_model_root, self.nlte.nlte_atom_data)

    def time_plasma_nlte_section_lu_config(self):
        self.nlte.tardis_config_verysimple_nlte["plasma"]["continuum_interaction"][
            "species"
        ] = [
            "He I",
        ]
        self.nlte.tardis_config_verysimple_nlte["plasma"]["nlte_ionization_species"] = ["H I"]
        self.nlte.tardis_config_verysimple_nlte["plasma"]["nlte_solver"] = "lu"
        config = Configuration.from_config_dict(self.nlte.tardis_config_verysimple_nlte)
        assemble_plasma(config, self.nlte.nlte_raw_model_lu, self.nlte.nlte_atom_data)

    def time_plasma_nlte_root_exc_section_config(self):
        self.nlte.tardis_config_verysimple_nlte["plasma"]["continuum_interaction"][
            "species"
        ] = [
            "He I",
        ]
        self.nlte.tardis_config_verysimple_nlte["plasma"]["nlte_excitation_species"] = ["H I"]
        self.nlte.tardis_config_verysimple_nlte["plasma"]["nlte_solver"] = "root"
        config = Configuration.from_config_dict(self.nlte.tardis_config_verysimple_nlte)
        assemble_plasma(config, self.nlte.nlte_raw_model_root, self.nlte.nlte_atom_data)

    def time_spectrum_section_config(self):
        self.tardis_config_verysimple["spectrum"]["start"] = Quantity("2500 angstrom")
        self.tardis_config_verysimple["spectrum"]["stop"] = Quantity("500 angstrom")
        Configuration.from_config_dict(
            self.tardis_config_verysimple, validate=True, config_dirname="test"
        )
