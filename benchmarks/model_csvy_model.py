"""
Basic TARDIS Benchmark.
"""
import numpy as np
import pandas as pd
from asv_runner.benchmarks.mark import parameterize, skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.io.configuration.config_reader import Configuration
from tardis.model import SimulationState


# @skip_benchmark
class BenchmarkModelCsvyModel(BenchmarkBase):
    """
    Class to benchmark the csvy model function.
    """

    def __init__(self):
        pass

    def model_config_fnames(self, model_name):
        csvy_config_file = f"{self.example_csvy_file_dir}/{model_name}_csvy.yml"
        old_config_file = f"{self.example_csvy_file_dir}/{model_name}_old_config.yml"
        return csvy_config_file, old_config_file

    @parameterize({"Model name": (
            "model_full", "branch85", "uniform", "powerlaw", "exponential", "radiative"
    )})
    def time_compare_models(self, model_name):
        csvy_config_file, old_config_file = self.model_config_fnames(model_name)
        tardis_config = Configuration.from_yaml(csvy_config_file)
        tardis_config_old = Configuration.from_yaml(old_config_file)
        csvy_simulation_state = SimulationState.from_csvy(
            tardis_config, atom_data=self.atomic_dataset
        )
        config_simulation_state = SimulationState.from_config(
            tardis_config_old, atom_data=self.atomic_dataset
        )
        csvy_simulation_state.get_properties().keys()
        config_simulation_state.get_properties().keys()

        assert (
                csvy_simulation_state.abundance.shape
                == config_simulation_state.abundance.shape
        )
        assert (
                csvy_simulation_state.composition.nuclide_mass_fraction.shape
                == config_simulation_state.composition.nuclide_mass_fraction.shape
        )
        assert (
                csvy_simulation_state.abundance.shape
                == config_simulation_state.abundance.shape
        )

    def time_dimensionality_after_update_v_inner_boundary(self):
        csvy_config_file = f"{self.example_csvy_file_dir}/radiative_csvy.yml"
        config = Configuration.from_yaml(csvy_config_file)
        csvy_model = SimulationState.from_csvy(config, atom_data=self.atomic_dataset)

        new_config = config
        new_config.model.v_inner_boundary = csvy_model.velocity[1]
        new_csvy_model = SimulationState.from_csvy(
            new_config, atom_data=self.atomic_dataset
        )

        assert new_csvy_model.no_of_raw_shells == csvy_model.no_of_raw_shells
        assert new_csvy_model.no_of_shells == csvy_model.no_of_shells - 1
        assert new_csvy_model.velocity.shape[0] == csvy_model.velocity.shape[0] - 1
        assert new_csvy_model.density.shape[0] == csvy_model.density.shape[0] - 1
        assert new_csvy_model.volume.shape[0] == csvy_model.volume.shape[0] - 1
        assert (
                new_csvy_model.t_radiative.shape[0]
                == csvy_model.t_radiative.shape[0] - 1
        )

    @property
    def csvy_model_test_abundances(self):
        csvypath = f"{self.example_csvy_file_dir}/csvy_model_to_test_abundances.yml"
        config = Configuration.from_yaml(csvypath)
        csvy_model_test_abundances = SimulationState.from_csvy(
            config, atom_data=self.atomic_dataset
        )
        return csvy_model_test_abundances

    @property
    def reference_input_dataframes(self):
        abundance_index = pd.Index([1, 2], name="atomic_number")
        reference_input_abundance = pd.DataFrame(
            [[0.0, 0.33, 0.3, 0.5, 0.4, 0.2], [0.98, 0.64, 0.6, 0.4, 0.55, 0.79]],
            index=abundance_index,
        )

        arrays = [[28], [56]]
        isotope_index = pd.MultiIndex.from_arrays(
            arrays, names=["atomic_number", "mass_number"]
        )
        reference_input_isotopes = pd.DataFrame(
            [[0.02, 0.03, 0.1, 0.1, 0.05, 0.01]],
            columns=np.arange(6),
            index=isotope_index,
        )
        return reference_input_abundance, reference_input_isotopes

    def time_read_csvy_abundances(self):
        (
            reference_input_abundance,
            reference_input_isotopes,
        ) = self.reference_input_dataframes

        composition = self.csvy_model_test_abundances.composition
        nuclide_mass_fraction = composition.nuclide_mass_fraction
        model_abundances = nuclide_mass_fraction[
            nuclide_mass_fraction.index.get_level_values(1) == -1
            ]

        reference_input_shape = reference_input_abundance.shape
        assert model_abundances.shape == reference_input_shape

    @property
    def reference_decayed_abundance(self):
        decay_index = pd.Index([1, 2, 26, 27, 28], name="atomic_number")
        reference_decayed_abundance = pd.DataFrame(
            [
                [0.0, 0.33, 0.3, 0.5, 0.4, 0.2],
                [0.98, 0.64, 0.6, 0.4, 0.55, 0.79],
                [
                    0.00013977843354947162,
                    0.00020966765032420787,
                    0.0006988921677473642,
                    0.0006988921677473642,
                    0.0003494460838736821,
                    6.988921677473581e-05,
                ],
                [
                    0.007188928223953217,
                    0.010783392335929825,
                    0.035944641119766085,
                    0.035944641119766085,
                    0.017972320559883043,
                    0.0035944641119766084,
                ],
                [
                    0.012671293342497312,
                    0.019006940013745966,
                    0.06335646671248656,
                    0.06335646671248656,
                    0.03167823335624328,
                    0.006335646671248656,
                ],
            ],
            index=decay_index,
        )
        return reference_decayed_abundance

    def time_csvy_model_decay(self):
        model_decayed_abundance_shape = self.csvy_model_test_abundances.abundance.shape
        reference_decayed_abundance_shape = self.reference_decayed_abundance.shape
        assert model_decayed_abundance_shape == reference_decayed_abundance_shape
