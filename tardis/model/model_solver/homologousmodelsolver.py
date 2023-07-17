from tardis.model.model_solver import ModelSolver
from tardis.model.density import HomologousDensityState


class HomologousModelStateSolver(ModelSolver):
    def evolve_geometry(self, new_time):
        geometry = self._model_state.geometry
        geometry.r_inner = geometry.v_inner * new_time
        geometry.r_outer = geometry.v_outer * new_time
        return geometry

    def evolve_density(self, new_time):
        density = self._model_state.composition.density
        new_density = density * (new_time / density.time) ** -3
        return HomologousDensityState(new_density, new_time)
