from collections import OrderedDict
from copy import deepcopy
from pathlib import Path

from benchmarks.util.base import Base
from tardis.io.atom_data import AtomData
from tardis.io.configuration.config_reader import Configuration
from tardis.io.util import yaml_load_file, YAMLLoader
from tardis.model import SimulationState


class NLTE:
    def __init__(self):
        self.__base = Base()

    @property
    def tardis_config_verysimple_nlte(self) -> OrderedDict:
        path: str = (
            "../../tardis/io/configuration/tests/data/tardis_configv1_nlte.yml"
        )
        filename: Path = self.__base.get_path(path)

        return yaml_load_file(
            filename,
            YAMLLoader,
        )

    @property
    def nlte_raw_model_root(self) -> SimulationState:
        return SimulationState.from_config(
            self.tardis_model_config_nlte_root, self.nlte_atom_data
        )

    @property
    def nlte_raw_model_lu(self) -> SimulationState:
        return SimulationState.from_config(
            self.tardis_model_config_nlte_lu, self.nlte_atom_data
        )

    @property
    def nlte_atom_data(self) -> AtomData:
        atomic_data = deepcopy(self.nlte_atomic_dataset)
        return atomic_data

    @property
    def nlte_atomic_dataset(self) -> AtomData:
        nlte_atomic_data = AtomData.from_hdf(self.nlte_atomic_data_fname)
        return nlte_atomic_data

    @property
    def nlte_atomic_data_fname(self) -> str:
        atomic_data_fname = (
            f"{self.__base.tardis_ref_path}/nlte_atom_data/TestNLTE_He_Ti.h5"
        )

        if not Path(atomic_data_fname).exists():
            atom_data_missing_str = (
                f"Atomic datafiles {atomic_data_fname} does not seem to exist"
            )
            raise Exception(atom_data_missing_str)

        return atomic_data_fname

    @property
    def tardis_model_config_nlte_root(self) -> Configuration:
        config = Configuration.from_yaml(
            f"{self.__base.example_configuration_dir}/tardis_configv1_nlte.yml"
        )
        config.plasma.nlte_solver = "root"
        return config

    @property
    def tardis_model_config_nlte_lu(self) -> Configuration:
        config = Configuration.from_yaml(
            f"{self.__base.example_configuration_dir}/tardis_configv1_nlte.yml"
        )
        config.plasma.nlte_solver = "lu"
        return config
