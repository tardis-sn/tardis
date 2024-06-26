from copy import deepcopy
from os.path import dirname, join, realpath
from pathlib import Path
from tempfile import mkstemp

import astropy.units as u
import numpy as np
import pandas as pd
from numba import njit

from benchmarks.util.nlte import NLTE
from tardis.io.atom_data import AtomData
from tardis.io.configuration import config_reader
from tardis.io.configuration.config_reader import Configuration
from tardis.io.util import YAMLLoader, yaml_load_file
from tardis.model import SimulationState
from tardis.simulation import Simulation
from tardis.tests.fixtures.atom_data import DEFAULT_ATOM_DATA_UUID
from tardis.tests.fixtures.regression_data import RegressionData
from tardis.transport.montecarlo import RPacket
from tardis.transport.montecarlo.numba_interface import (
    opacity_state_initialize,
)
from tardis.transport.montecarlo.packet_collections import (
    VPacketCollection,
)


class BenchmarkBase:
    # It allows 10 minutes of runtime for each benchmark and includes
    # the total time for all the repetitions for each benchmark.
    timeout = 600

    def __init__(self):
        self.nlte = NLTE()

    @staticmethod
    def get_relative_path(partial_path: str):
        path = dirname(realpath(__file__))
        targets = Path(partial_path).parts

        for target in targets:
            path = join(path, target)

        return path

    def get_absolute_path(self, partial_path):
        partial_path = "../" + partial_path

        return self.get_relative_path(partial_path)

    @property
    def tardis_config_verysimple(self):
        filename = self.get_absolute_path(
            "tardis/io/configuration/tests/data/tardis_configv1_verysimple.yml"
        )
        return yaml_load_file(
            filename,
            YAMLLoader,
        )

    @property
    def tardis_ref_path(self):
        # TODO: This route is fixed but needs to get from the arguments given in the command line.
        #       /app/tardis-refdata
        ref_data_path = Path(
            Path(__file__).parent.parent,
            "tardis-refdata",
        ).resolve()
        return ref_data_path

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
        atomic_data_fname = (
            f"{self.tardis_ref_path}/atom_data/kurucz_cd23_chianti_H_He.h5"
        )

        if not Path(atomic_data_fname).exists():
            atom_data_missing_str = (
                f"{atomic_data_fname} atomic datafiles "
                f"does not seem to exist"
            )
            raise Exception(atom_data_missing_str)

        return atomic_data_fname

    @property
    def example_configuration_dir(self):
        return self.get_absolute_path("tardis/io/configuration/tests/data")

    @property
    def hdf_file_path(self):
        # TODO: Delete this files after ASV runs the benchmarks.
        #       ASV create a temporal directory in runtime per test: `tmpiuxngvlv`.
        #       The ASV and ASV_Runner, not has some way to get this temporal directory.
        #       The idea is use this temporal folders to storage this created temporal file.
        _, path = mkstemp("-tardis-benchmark-hdf_buffer-test.hdf")
        return path

    def create_temporal_file(self, suffix=None):
        # TODO: Delete this files after ASV runs the benchmarks.
        #       ASV create a temporal directory in runtime per test: `tmpiuxngvlv`.
        #       The ASV and ASV_Runner, not has some way to get this temporal directory.
        #       The idea is use this temporal folders to storage this created temporal file.
        suffix_str = "" if suffix is None else f"-{suffix}"
        _, path = mkstemp(suffix_str)
        return path

    @property
    def gamma_ray_simulation_state(self):
        self.gamma_ray_config.model.structure.velocity.start = 1.0 * u.km / u.s
        self.gamma_ray_config.model.structure.density.rho_0 = (
            5.0e2 * u.g / u.cm**3
        )
        self.gamma_ray_config.supernova.time_explosion = 150 * u.d

        return SimulationState.from_config(
            self.gamma_ray_config, atom_data=self.atomic_dataset
        )

    @property
    def gamma_ray_config(self):
        yml_path = f"{self.example_configuration_dir}/tardis_configv1_density_exponential_nebular_multi_isotope.yml"

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
        sim = Simulation.from_config(
            self.config_verysimple, atom_data=atomic_data
        )
        sim.iterate(4000)
        return sim

    @property
    def config_verysimple(self):
        return Configuration.from_yaml(
            f"{self.example_configuration_dir}/tardis_configv1_verysimple.yml"
        )

    class CustomPyTestRequest:
        def __init__(
            self,
            tardis_regression_data_path: str,
            node_name: str,
            node_module_name: str,
            regression_data_dir: str,
        ):
            self.tardis_regression_data_path = tardis_regression_data_path
            self.node_name = node_name
            self.node_module_name = node_module_name
            self.regression_data_dir = regression_data_dir

        @property
        def config(self):
            class SubClass:
                @staticmethod
                def getoption(option):
                    if option == "--tardis-regression-data":
                        return self.tardis_regression_data_path
                    return None

            return SubClass()

        @property
        def node(self):
            class SubClass:
                def __init__(self, parent):
                    self.parent = parent

                @property
                def name(self):
                    return self.parent.node_name

                @property
                def module(self):
                    class SubSubClass:
                        def __init__(self, parent):
                            self.parent = parent

                        @property
                        def __name__(self):
                            return self.parent.node_module_name

                    return SubSubClass(self.parent)

            return SubClass(self)

        @property
        def cls(self):
            return None

        @property
        def relative_regression_data_dir(self):
            return self.regression_data_dir

    @staticmethod
    def regression_data(request: CustomPyTestRequest):
        return RegressionData(request)

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
        return (
            self.nb_simulation_verysimple.transport.transport_state.packet_collection
        )

    @property
    def nb_simulation_verysimple(self):
        atomic_data = deepcopy(self.atomic_dataset)
        sim = Simulation.from_config(
            self.config_verysimple, atom_data=atomic_data
        )
        sim.iterate(10)
        return sim

    @property
    def verysimple_time_explosion(self):
        model = self.nb_simulation_verysimple.simulation_state
        return model.time_explosion.cgs.value

    @property
    def verysimple_opacity_state(self):
        return opacity_state_initialize(
            self.nb_simulation_verysimple.plasma,
            line_interaction_type="macroatom",
            disable_line_scattering=self.nb_simulation_verysimple.transport.montecarlo_configuration.DISABLE_LINE_SCATTERING,
            continuum_processes_enabled=self.nb_simulation_verysimple.transport.montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED,
        )

    @property
    def verysimple_enable_full_relativity(self):
        return self.nb_simulation_verysimple.transport.enable_full_relativity

    @property
    def verysimple_disable_line_scattering(self):
        return (
            self.nb_simulation_verysimple.transport.montecarlo_configuration.DISABLE_LINE_SCATTERING
        )

    @property
    def verysimple_continuum_processes_enabled(self):
        return (
            self.nb_simulation_verysimple.transport.montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED
        )

    @property
    def verysimple_tau_russian(self):
        return (
            self.nb_simulation_verysimple.transport.montecarlo_configuration.VPACKET_TAU_RUSSIAN
        )

    @property
    def verysimple_survival_probability(self):
        return (
            self.nb_simulation_verysimple.transport.montecarlo_configuration.SURVIVAL_PROBABILITY
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
        return (
            self.nb_simulation_verysimple.simulation_state.geometry.to_numba()
        )

    @property
    def simulation_verysimple_vpacket_tracking(self):
        atomic_data = deepcopy(self.atomic_dataset)
        sim = Simulation.from_config(
            self.config_verysimple,
            atom_data=atomic_data,
            virtual_packet_logging=True,
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
        with pd.HDFStore(
            f"{self.tardis_ref_path}/unit_test_data.h5", mode=mode
        ) as store:
            yield store
