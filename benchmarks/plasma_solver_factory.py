"""
Basic TARDIS Benchmark.
"""

from copy import deepcopy

from benchmarks.benchmark_base import BenchmarkBase
from tardis.plasma.assembly.base import PlasmaSolverFactory
from tardis.plasma.radiation_field import DilutePlanckianRadiationField


class BenchmarkPlasmaSolverFactory(BenchmarkBase):
    """
    Class to benchmark the PlasmaSolverFactory.
    """

    repeat = 2

    def setup(self):
        """Set up simulation inputs for benchmarking."""
        self.sim = self.nb_simulation_verysimple
        self.config = self.config_verysimple

        sim_state = self.sim.simulation_state
        self.atomic_numbers = sim_state.abundance.index

        atom_data = deepcopy(self.atomic_dataset)
        self.prepared_factory = PlasmaSolverFactory(atom_data, self.config)
        self.prepared_factory.prepare_factory(
            self.atomic_numbers,
            "tardis.plasma.properties.legacy_property_collections",
            self.config,
        )
        self.number_densities = sim_state.calculate_elemental_number_density(
            atom_data.atom_data.mass.copy()
        )
        self.dilute_planckian_radiation_field = DilutePlanckianRadiationField(
            sim_state.t_radiative, sim_state.dilution_factor
        )
        self.time_explosion = sim_state.time_explosion
        self.electron_densities = sim_state._electron_densities

    def time_prepare_factory(self):
        """Benchmark factory creation and module resolution."""
        atom_data = deepcopy(self.atomic_dataset)
        factory = PlasmaSolverFactory(atom_data, self.config)
        factory.prepare_factory(
            self.atomic_numbers,
            "tardis.plasma.properties.legacy_property_collections",
            self.config,
        )

    def time_assemble(self):
        """Benchmark plasma assembly from a prepared factory."""
        self.prepared_factory.assemble(
            self.number_densities,
            self.dilute_planckian_radiation_field,
            self.time_explosion,
            self.electron_densities,
        )
