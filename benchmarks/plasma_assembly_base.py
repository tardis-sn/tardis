

import functools

from benchmarks.benchmark_base import BenchmarkBase
from tardis.plasma.assembly import PlasmaSolverFactory
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


class BenchmarkPlasmaAssemblyBase(BenchmarkBase):
    repeat = 1

    @functools.cache
    def setup(self):
        self.PlasmaSolverFactory = PlasmaSolverFactory(
            self.atomic_dataset,
            self.config_verysimple,
        )
        self.atomic_numbers = self.simulation_state.abundance.index
        self.dilute_planckian_radiation_field = DilutePlanckianRadiationField(
            self.simulation_state.t_radiative, self.simulation_state.dilution_factor
        )
        

    def time_prepare_factory(self):
        self.PlasmaSolverFactory.prepare_factory( property_collections="tardis.plasma.properties.legacy_property_collections", selected_atomic_numbers=self.nb_simulation_verysimple.simulation_state.abundance.index)

    def time_setup_helium_treatment(self):
        self.PlasmaSolverFactory.setup_helium_treatment()

    def time_setup_analytical_approximations(self):
        self.PlasmaSolverFactory.setup_analytical_approximations()

    def time_setup_analytical_approximations(self):
        self.PlasmaSolverFactory.setup_analytical_approximations()

    def time_assemble(self):
        self.PlasmaSolverFactory.prepare_factory(
            self.atomic_numbers,
            "tardis.plasma.properties.legacy_property_collections",
            self.config_verysimple,
        )
        self.PlasmaSolverFactory.assemble(
            self.simulation_state.elemental_number_density,
            self.dilute_planckian_radiation_field,
            self.simulation_state.time_explosion,
            self.simulation_state._electron_densities,
        )
    
    def time_setup_continuum_interactions(self):
        self.PlasmaSolverFactory.setup_continuum_interactions()

    def time_initialize_j_blues(self):
        self.PlasmaSolverFactory.initialize_j_blues(
            self.dilute_planckian_radiation_field,
            self.PlasmaSolverFactory.atom_data.lines
        )
