"""
Basic TARDIS Benchmark.
"""

from copy import deepcopy

from astropy import units as u

from benchmarks.benchmark_base import BenchmarkBase
from tardis.plasma.assembly.base import PlasmaSolverFactory


class BenchmarkPlasmaSolverFactory(BenchmarkBase):
    """
    Benchmarks for the `PlasmaSolverFactory` workflow.
    """

    repeat = 2

    def setup(self):
        self.config = self.config_verysimple
        self.nb_simulation = self.nb_simulation_verysimple

        self.electron_densities = (
            self.nb_simulation.plasma.electron_densities.values * (u.cm ** -3)
        )
        self.selected_atomic_numbers = (
            self.nb_simulation.plasma.atomic_data.selected_atomic_numbers
        )
        self.property_collections = (
            "tardis.plasma.properties.legacy_property_collections"
        )

        self.number_density = self.nb_simulation.plasma.number_density * (u.cm ** -3)
        self.dilute_planckian_radiation_field = (
            self.nb_simulation.plasma.dilute_planckian_radiation_field
        )
        self.time_explosion = self.nb_simulation.plasma.time_explosion

        atom_data = deepcopy(self.atomic_dataset)
        self.solver = PlasmaSolverFactory(atom_data, self.config)
        self.solver.prepare_factory(
            self.selected_atomic_numbers,
            self.property_collections,
            self.config,
        )

    def time_prepare_factory(self):
        atom_data = deepcopy(self.atomic_dataset)
        solver = PlasmaSolverFactory(atom_data, self.config)
        solver.prepare_factory(
            self.selected_atomic_numbers,
            self.property_collections,
            self.config,
        )

    def time_assemble(self):
        self.solver.assemble(
            self.number_density,
            self.dilute_planckian_radiation_field,
            self.time_explosion,
            self.electron_densities,
        )
