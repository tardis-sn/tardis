import re
from copy import deepcopy
from os.path import dirname, realpath, join
from pathlib import Path
from tempfile import mkstemp

import astropy.units as u
import numpy as np
import pandas as pd
from numba import njit

from tardis.io.atom_data import AtomData
from tardis.io.configuration import config_reader
from tardis.io.configuration.config_reader import Configuration
from tardis.io.util import yaml_load_file, YAMLLoader, HDFWriterMixin
from tardis.model import SimulationState
from tardis.montecarlo import NumbaModel, opacity_state_initialize
from tardis.montecarlo.montecarlo_numba import RPacket
from tardis.montecarlo.montecarlo_numba.packet_collections import (
    VPacketCollection,
)
from tardis.simulation import Simulation
from tardis.tests.fixtures.atom_data import DEFAULT_ATOM_DATA_UUID


class BenchmarkBase:
    timeout = 369

    @staticmethod
    def split_path(path):
        return path.split('/')

    def get_relative_path(self, partial_path: str):
        path = dirname(realpath(__file__))
        targets = self.split_path(partial_path)

        for target in targets:
            path = join(path, target)

        return path

    def get_absolute_path(self, partial_path):
        partial_path = "../" + partial_path

        return self.get_relative_path(partial_path)

    @property
    def tardis_config_verysimple(self):
        filename = self.get_absolute_path("tardis/io/configuration/tests/data/tardis_configv1_verysimple.yml")
        return yaml_load_file(
            filename,
            YAMLLoader,
        )

    @property
    def tardis_config_verysimple_nlte(self):
        filename = self.get_absolute_path("tardis/io/configuration/tests/data/tardis_configv1_nlte.yml")
        return yaml_load_file(
            filename,
            YAMLLoader,
        )

    @property
    def nlte_raw_model_root(self):
        return SimulationState.from_config(
            self.tardis_model_config_nlte_root, self.nlte_atom_data
        )

    @property
    def nlte_raw_model_lu(self):
        return SimulationState.from_config(
            self.tardis_model_config_nlte_lu, self.nlte_atom_data
        )

    @property
    def nlte_atom_data(self):
        atomic_data = deepcopy(self.nlte_atomic_dataset)
        return atomic_data

    @property
    def nlte_atomic_dataset(self):
        nlte_atomic_data = AtomData.from_hdf(self.nlte_atomic_data_fname)
        return nlte_atomic_data

    @property
    def nlte_atomic_data_fname(self):
        atomic_data_fname = f"{self.tardis_ref_path}/nlte_atom_data/TestNLTE_He_Ti.h5"
        atom_data_missing_str = f"{atomic_data_fname} atomic datafiles " f"does not seem to exist"

        if not Path(atomic_data_fname).exists():
            raise Exception(atom_data_missing_str)

        return atomic_data_fname

    @property
    def tardis_ref_path(self):
        # TODO: This route is fixed but needs to get from the arguments given in the command line.
        #       /app/tardis-refdata
        return '/app/tardis-refdata'

    @property
    def atomic_dataset(self) -> AtomData:
        atomic_data = AtomData.from_hdf(self.atomic_data_fname)

        if atomic_data.md5 != DEFAULT_ATOM_DATA_UUID:
            message = f'Need default Kurucz atomic dataset (md5="{DEFAULT_ATOM_DATA_UUID}")'
            raise Exception(message)
        else:
            return atomic_data

    @property
    def atomic_data_fname(self):
        atomic_data_fname = f"{self.tardis_ref_path}/atom_data/kurucz_cd23_chianti_H_He.h5"

        if not Path(atomic_data_fname).exists():
            atom_data_missing_str = f"{atomic_data_fname} atomic datafiles " f"does not seem to exist"
            raise Exception(atom_data_missing_str)

        return atomic_data_fname

    @property
    def tardis_model_config_nlte_root(self):
        config = Configuration.from_yaml(
            f"{self.example_configuration_dir}/tardis_configv1_nlte.yml"
        )
        config.plasma.nlte_solver = "root"
        return config

    @property
    def tardis_model_config_nlte_lu(self):
        config = Configuration.from_yaml(
            f"{self.example_configuration_dir}/tardis_configv1_nlte.yml"
        )
        config.plasma.nlte_solver = "lu"
        return config

    @property
    def example_configuration_dir(self):
        return self.get_absolute_path("tardis/io/configuration/tests/data")

    @property
    def hdf_file_path(self):
        # TODO: Delete this files after ASV runs the benchmarks.
        #       ASV create a temporal directory in runtime per test: `tmpiuxngvlv`.
        #       The ASV and ASV_Runner, not has some way to get this temporal directory.
        #       The idea is use this temporal folders to storage this created temporal file.
        _, path = mkstemp('-tardis-benchmark-hdf_buffer-test.hdf')
        return path

    def create_temporal_file(self, suffix=None):
        # TODO: Delete this files after ASV runs the benchmarks.
        #       ASV create a temporal directory in runtime per test: `tmpiuxngvlv`.
        #       The ASV and ASV_Runner, not has some way to get this temporal directory.
        #       The idea is use this temporal folders to storage this created temporal file.
        suffix_str = '' if suffix is None else f'-{suffix}'
        _, path = mkstemp(suffix_str)
        return path

    @property
    def gamma_ray_simulation_state(self):
        self.gamma_ray_config.model.structure.velocity.start = 1.0 * u.km / u.s
        self.gamma_ray_config.model.structure.density.rho_0 = 5.0e2 * u.g / u.cm ** 3
        self.gamma_ray_config.supernova.time_explosion = 150 * u.d

        return SimulationState.from_config(
            self.gamma_ray_config, atom_data=self.atomic_dataset
        )

    @property
    def gamma_ray_config(self):
        yml_path = (
            f"{self.example_configuration_dir}/tardis_configv1_density_exponential_nebular_multi_isotope.yml"
        )

        return config_reader.Configuration.from_yaml(yml_path)

    @property
    def example_model_file_dir(self):
        return self.get_absolute_path("tardis/io/model/readers/tests/data")

    @property
    def kurucz_atomic_data(self) -> AtomData:
        return deepcopy(self.atomic_dataset)

    @property
    def example_csvy_file_dir(self):
        return self.get_absolute_path("tardis/model/tests/data/")

    @property
    def simulation_verysimple(self):
        atomic_data = deepcopy(self.atomic_dataset)
        sim = Simulation.from_config(self.config_verysimple, atom_data=atomic_data)
        sim.iterate(4000)
        return sim

    @property
    def config_verysimple(self):
        return Configuration.from_yaml(
            f"{self.example_configuration_dir}/tardis_configv1_verysimple.yml"
        )

    class RegressionData:
        def __init__(self) -> None:
            # TODO: This route is fixed but needs to get from the arguments given in the command line.
            #       /app/tardis-regression-data
            self.regression_data_path = "/app/tardis-regression-data"
            # TODO: Parameter `--generate-reference` is set in the command line.
            #       It is fixed but needs to get from the arguments.
            self.enable_generate_reference = False
            self.fname = f"{self.fname_prefix}.UNKNOWN_FORMAT"

        @property
        def module_name(self):
            return self.__name__

        @property
        def test_name(self):
            return self.name

        @property
        def fname_prefix(self):
            double_under = re.compile(r"[:\[\]{}]")
            no_space = re.compile(r'[,"\']')  # quotes and commas

            name = double_under.sub("__", self.test_name)
            name = no_space.sub("", name)
            return name

        @property
        def relative_regression_data_dir(self):
            relative_data_dir = Path(self.module_name.replace(".", "/"))
            if self.cls is not None:
                relative_data_dir /= HDFWriterMixin.convert_to_snake_case(
                    self.cls.__name__
                )
            return relative_data_dir

        @property
        def absolute_regression_data_dir(self):
            return f"{self.regression_data_path}/{self.relative_regression_data_dir}"

        @property
        def fpath(self):
            return f"{self.absolute_regression_data_dir}/{self.fname}"

        def sync_dataframe(self, data, key="data"):
            """
            Synchronizes the dataframe with the regression data.

            Parameters
            ----------
            data : DataFrame
                The dataframe to be synchronized.
            key : str, optional
                The key to use for storing the dataframe in the regression data file. Defaults to "data".

            Returns
            -------
            DataFrame or None
                The synchronized dataframe if `enable_generate_reference` is `False`, otherwise `None`.
            """
            self.fname = f"{self.fname_prefix}.h5"
            if self.enable_generate_reference:
                Path(self.fpath).parent.mkdir(parents=True, exist_ok=True)
                data.to_hdf(
                    self.fpath,
                    key=key,
                )
                raise Exception("Skipping test to generate reference data")
            else:
                return pd.read_hdf(self.fpath, key=key)

        def sync_ndarray(self, data):
            """
            Synchronizes the ndarray with the regression data.

            Parameters
            ----------
            data : ndarray
                The ndarray to be synchronized.

            Returns
            -------
            ndarray or None
                The synchronized ndarray if `enable_generate_reference` is `False`, otherwise `None`.
            """
            self.fname = f"{self.fname_prefix}.npy"
            if self.enable_generate_reference:
                Path(self.fpath).parent.mkdir(parents=True, exist_ok=True)
                Path(self.fpath).parent.mkdir(parents=True, exist_ok=True)
                np.save(self.fpath, data)
                raise Exception("Skipping test to generate reference data")
            else:
                return np.load(self.fpath)

        def sync_str(self, data):
            """
            Synchronizes the string with the regression data.

            Parameters
            ----------
            data : str
                The string to be synchronized.

            Returns
            -------
            str or None
                The synchronized string if `enable_generate_reference` is `False`, otherwise `None`.
            """
            self.fname = f"{self.fname_prefix}.txt"
            if self.enable_generate_reference:
                Path(self.fpath).parent.mkdir(parents=True, exist_ok=True)
                with Path(self.fpath).open("w") as fh:
                    fh.write(data)
                raise Exception(
                    f"Skipping test to generate regression_data {self.fpath} data"
                )
            else:
                with Path(self.fpath).open("r") as fh:
                    return fh.read()

        def sync_hdf_store(self, tardis_module, update_fname=True):
            """
            Synchronizes the HDF store with the regression data.

            Parameters
            ----------
            tardis_module : object
                The module to be synchronized.
            update_fname : bool, optional
                Whether to update the file name. Defaults to True.

            Returns
            -------
            HDFStore or None
                The synchronized HDF store if `enable_generate_reference` is `False`, otherwise `None`.
            """
            if update_fname:
                self.fname = f"{self.fname_prefix}.h5"
            if self.enable_generate_reference:
                Path(self.fpath).parent.mkdir(parents=True, exist_ok=True)
                with pd.HDFStore(self.fpath, mode="w") as store:
                    tardis_module.to_hdf(store, overwrite=True)
                raise Exception(
                    f"Skipping test to generate regression_data {self.fpath} data"
                )
            else:
                return pd.HDFStore(self.fpath, mode="r")

    @property
    def regression_data(self):
        return self.RegressionData()

    @property
    def packet(self):
        return RPacket(
            r=7.5e14,
            nu=self.verysimple_packet_collection.initial_nus[0],
            mu=self.verysimple_packet_collection.initial_mus[0],
            energy=self.verysimple_packet_collection.initial_energies[0],
            seed=1963,
            index=0,
        )

    @property
    def verysimple_packet_collection(self):
        return self.nb_simulation_verysimple.transport.transport_state.packet_collection

    @property
    def nb_simulation_verysimple(self):
        atomic_data = deepcopy(self.atomic_dataset)
        sim = Simulation.from_config(self.config_verysimple, atom_data=atomic_data)
        sim.iterate(10)
        return sim

    @property
    def verysimple_numba_model(self):
        model = self.nb_simulation_verysimple.simulation_state
        return NumbaModel(
            model.time_explosion.to("s").value,
        )

    @property
    def verysimple_opacity_state(self):
        return opacity_state_initialize(
            self.nb_simulation_verysimple.plasma, line_interaction_type="macroatom"
        )

    @property
    def static_packet(self):
        return RPacket(
            r=7.5e14,
            nu=0.4,
            mu=0.3,
            energy=0.9,
            seed=1963,
            index=0,
        )

    @property
    def set_seed_fixture(self):
        def set_seed(value):
            np.random.seed(value)

        return njit(set_seed)

    @property
    def verysimple_3vpacket_collection(self):
        spectrum_frequency = (
            self.nb_simulation_verysimple.transport.spectrum_frequency.value
        )
        return VPacketCollection(
            source_rpacket_index=0,
            spectrum_frequency=spectrum_frequency,
            number_of_vpackets=3,
            v_packet_spawn_start_frequency=0,
            v_packet_spawn_end_frequency=np.inf,
            temporary_v_packet_bins=0,
        )

    @property
    def verysimple_numba_radial_1d_geometry(self):
        return self.nb_simulation_verysimple.simulation_state.geometry.to_numba()

    @property
    def simulation_verysimple_vpacket_tracking(self):
        atomic_data = deepcopy(self.atomic_dataset)
        sim = Simulation.from_config(
            self.config_verysimple, atom_data=atomic_data, virtual_packet_logging=True
        )
        sim.last_no_of_packets = 4000
        sim.run_final()
        return sim

    @property
    def generate_reference(self):
        # TODO: Investigate how to get the `--generate-reference` parameter passed in the command line.
        #       `request.config.getoption("--generate-reference")`
        option = None
        if option is None:
            return False
        else:
            return option

    @property
    def tardis_ref_data(self):
        # TODO: This function is not working in the benchmarks.
        if self.generate_reference:
            mode = "w"
        else:
            mode = "r"
        with pd.HDFStore(f"{self.tardis_ref_path}/unit_test_data.h5", mode=mode) as store:
            yield store
