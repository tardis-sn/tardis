import functools
from os import environ as env
from copy import deepcopy
from pathlib import Path

import numpy as np

from tardis import run_tardis
from tardis.io.atom_data import AtomData
from tardis.io.configuration.config_reader import Configuration
from tardis.io.util import YAMLLoader, yaml_load_file
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.simulation import Simulation
from tardis.tests.fixtures.atom_data import DEFAULT_ATOM_DATA_MD5
from tardis.transport.montecarlo import RPacket, packet_trackers
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.estimators import radfield_mc_estimators
from tardis.opacities.opacity_state_numba import opacity_state_numba_initialize
from tardis.transport.montecarlo.packet_collections import VPacketCollection


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
    def tardis_config_verysimple(self):
        filename = self.get_absolute_path(
            "tardis/io/configuration/tests/data/tardis_configv1_verysimple.yml"
        )
        return yaml_load_file(
            filename,
            YAMLLoader,
        )

    @functools.cached_property
    def config_rpacket_tracking(self):
        config = Configuration.from_yaml(
            f"{self.example_configuration_dir}/tardis_configv1_verysimple.yml"
        )
        config.montecarlo.tracking.track_rpacket = True
        return config

    @functools.cached_property
    def tardis_ref_path(self):
        ref_data_path = Path(
            Path(__file__).parent.parent, env.get("TARDIS_REF_PATH")
        ).resolve()
        return ref_data_path

    @functools.cached_property
    def atomic_dataset(self) -> AtomData:
        atomic_data = AtomData.from_hdf(self.atomic_data_fname)

        if atomic_data.md5 != DEFAULT_ATOM_DATA_MD5:
            message = f'Need default Kurucz atomic dataset (md5="{DEFAULT_ATOM_DATA_MD5}")'
            raise Exception(message)
        else:
            return atomic_data

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
    def simulation_verysimple(self):
        atomic_data = deepcopy(self.atomic_dataset)
        sim = Simulation.from_config(
            self.config_verysimple, atom_data=atomic_data
        )
        sim.iterate(4000)
        return sim

    @functools.cached_property
    def config_verysimple(self):
        return Configuration.from_yaml(
            f"{self.example_configuration_dir}/tardis_configv1_verysimple.yml"
        )

    @functools.cached_property
    def packet(self):
        return RPacket(
            r=7.5e14,
            nu=self.verysimple_packet_collection.initial_nus[0],
            mu=self.verysimple_packet_collection.initial_mus[0],
            energy=self.verysimple_packet_collection.initial_energies[0],
            seed=1963,
            index=0,
        )

    @functools.cached_property
    def verysimple_packet_collection(self):
        return self.nb_simulation_verysimple.transport.transport_state.packet_collection

    @functools.cached_property
    def nb_simulation_verysimple(self):
        atomic_data = deepcopy(self.atomic_dataset)
        sim = Simulation.from_config(
            self.config_verysimple, atom_data=atomic_data
        )
        sim.iterate(5)
        return sim

    @functools.cached_property
    def verysimple_time_explosion(self):
        model = self.nb_simulation_verysimple.simulation_state
        return model.time_explosion.cgs.value

    @functools.cached_property
    def verysimple_opacity_state(self):
        return opacity_state_numba_initialize(
            self.nb_simulation_verysimple.plasma,
            line_interaction_type="macroatom",
            disable_line_scattering=self.nb_simulation_verysimple.transport.montecarlo_configuration.DISABLE_LINE_SCATTERING,
        )

    @functools.cached_property
    def verysimple_enable_full_relativity(self):
        return self.nb_simulation_verysimple.transport.enable_full_relativity

    @functools.cached_property
    def verysimple_tau_russian(self):
        return self.nb_simulation_verysimple.transport.montecarlo_configuration.VPACKET_TAU_RUSSIAN

    @functools.cached_property
    def verysimple_survival_probability(self):
        return self.nb_simulation_verysimple.transport.montecarlo_configuration.SURVIVAL_PROBABILITY

    @functools.cached_property
    def static_packet(self):
        return RPacket(
            r=7.5e14,
            nu=0.4,
            mu=0.3,
            energy=0.9,
            seed=1963,
            index=0,
        )

    @functools.cached_property
    def verysimple_3vpacket_collection(self):
        spectrum_frequency_grid = self.nb_simulation_verysimple.transport.spectrum_frequency_grid.value
        return VPacketCollection(
            source_rpacket_index=0,
            spectrum_frequency_grid=spectrum_frequency_grid,
            number_of_vpackets=3,
            v_packet_spawn_start_frequency=0,
            v_packet_spawn_end_frequency=np.inf,
            temporary_v_packet_bins=0,
        )

    @functools.cached_property
    def verysimple_numba_radial_1d_geometry(self):
        return (
            self.nb_simulation_verysimple.simulation_state.geometry.to_numba()
        )

    @functools.cached_property
    def montecarlo_configuration(self):
        return MonteCarloConfiguration()

    @functools.cached_property
    def rpacket_tracker(self):
        # Do not use RPacketTracker or RPacketLastInteraction directly
        # Use it by importing packet_trackers
        # functions with name track_* function is used by ASV
        return packet_trackers.RPacketLastInteractionTracker()

    @functools.cached_property
    def transport_state(self):
        return self.nb_simulation_verysimple.transport.transport_state

    @functools.cached_property
    def simulation_rpacket_tracking_enabled(self):
        config_verysimple = self.config_verysimple
        config_verysimple.montecarlo.iterations = 3
        config_verysimple.montecarlo.no_of_packets = 4000
        config_verysimple.montecarlo.last_no_of_packets = -1
        config_verysimple.spectrum.virtual.virtual_packet_logging = True
        config_verysimple.montecarlo.no_of_virtual_packets = 1
        config_verysimple.montecarlo.tracking.track_rpacket = True
        config_verysimple.spectrum.num = 2000
        atomic_data = deepcopy(self.atomic_dataset)
        sim = run_tardis(
            config_verysimple,
            atom_data=atomic_data,
            show_convergence_plots=False,
        )
        return sim

    @functools.cached_property
    def geometry(self):
        return NumbaRadial1DGeometry(
            r_inner=np.array([6.912e14, 8.64e14], dtype=np.float64),
            r_outer=np.array([8.64e14, 1.0368e15], dtype=np.float64),
            v_inner=np.array([-1, -1], dtype=np.float64),
            v_outer=np.array([-1, -1], dtype=np.float64),
        )

    @functools.cached_property
    def estimators(self):
        return radfield_mc_estimators.RadiationFieldMCEstimators(
            j_estimator=np.array([0.0, 0.0], dtype=np.float64),
            nu_bar_estimator=np.array([0.0, 0.0], dtype=np.float64),
            j_blue_estimator=np.array(
                [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=np.float64
            ),
            Edotlu_estimator=np.array(
                [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]], dtype=np.float64
            ),
            photo_ion_estimator=np.empty((0, 0), dtype=np.float64),
            stim_recomb_estimator=np.empty((0, 0), dtype=np.float64),
            bf_heating_estimator=np.empty((0, 0), dtype=np.float64),
            stim_recomb_cooling_estimator=np.empty((0, 0), dtype=np.float64),
            photo_ion_estimator_statistics=np.empty((0, 0), dtype=np.int64),
        )

    @functools.cached_property
    def rpacket_tracker_list(self):
        no_of_packets = len(self.transport_state.packet_collection.initial_nus)
        return packet_trackers.generate_rpacket_last_interaction_tracker_list(
            no_of_packets
        )
