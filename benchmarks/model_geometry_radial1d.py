"""
Basic TARDIS Benchmark.
"""
import numpy as np
from astropy import units as u
from asv_runner.benchmarks.mark import skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.model.geometry.radial1d import HomologousRadial1DGeometry


# @skip_benchmark
class BenchmarkModelGeometryRadial1d(BenchmarkBase):
    """
    Class to benchmark the radial 1D function.
    """

    def __init__(self):
        pass

    @property
    def homologous_radial1d_geometry(self):
        velocity = np.arange(8000, 21000, 1000) * u.km / u.s
        v_inner = velocity[:-1]
        v_outer = velocity[1:]
        time_explosion = 5 * u.day
        geometry = HomologousRadial1DGeometry(
            v_inner, v_outer, v_inner[0], v_outer[-1], time_explosion
        )
        return geometry

    def time_vb_indices(self):
        # Testing if the indices returned are correct when inner and outer
        # boundary are on the innermost and outermost shell

        homologous_radial1d_geometry = self.homologous_radial1d_geometry

        homologous_radial1d_geometry.v_inner_boundary = (
            homologous_radial1d_geometry.v_inner[0]
        )
        homologous_radial1d_geometry.v_outer_boundary = (
            homologous_radial1d_geometry.v_outer[-1]
        )
        assert homologous_radial1d_geometry.v_inner_boundary_index == 0
        assert homologous_radial1d_geometry.v_outer_boundary_index == len(
            homologous_radial1d_geometry.v_inner
        )
        vib_index = homologous_radial1d_geometry.v_inner_boundary_index
        vob_index = homologous_radial1d_geometry.v_outer_boundary_index
        assert np.all(
            homologous_radial1d_geometry.v_inner[vib_index:vob_index]
            == homologous_radial1d_geometry.v_inner
        )
        EPSILON_VELOCITY_SHIFT = 1 * u.km / u.s
        # pivoting around the inner boundary of the simulation

        homologous_radial1d_geometry.v_inner_boundary = (
                homologous_radial1d_geometry.v_inner[0] + EPSILON_VELOCITY_SHIFT
        )
        assert homologous_radial1d_geometry.v_inner_boundary_index == 0
        homologous_radial1d_geometry.v_inner_boundary = (
                homologous_radial1d_geometry.v_inner[0] - EPSILON_VELOCITY_SHIFT
        )
        assert homologous_radial1d_geometry.v_inner_boundary_index == 0

        # pivoting around the first shell boundary of the simulation
        homologous_radial1d_geometry.v_inner_boundary = (
                homologous_radial1d_geometry.v_inner[1] - EPSILON_VELOCITY_SHIFT
        )
        assert homologous_radial1d_geometry.v_inner_boundary_index == 0
        homologous_radial1d_geometry.v_inner_boundary = (
                homologous_radial1d_geometry.v_inner[1] + EPSILON_VELOCITY_SHIFT
        )
        assert homologous_radial1d_geometry.v_inner_boundary_index == 1

        # pivoting around the outer boundary of the simulation
        homologous_radial1d_geometry.v_outer_boundary = (
                homologous_radial1d_geometry.v_outer[-1] + EPSILON_VELOCITY_SHIFT
        )
        assert homologous_radial1d_geometry.v_outer_boundary_index == 12
        homologous_radial1d_geometry.v_outer_boundary = (
                homologous_radial1d_geometry.v_outer[-1] - EPSILON_VELOCITY_SHIFT
        )
        assert homologous_radial1d_geometry.v_outer_boundary_index == 12

        # pivoting around the second to outer boundary of the simulation
        homologous_radial1d_geometry.v_outer_boundary = (
                homologous_radial1d_geometry.v_outer[-2] + EPSILON_VELOCITY_SHIFT
        )
        assert homologous_radial1d_geometry.v_outer_boundary_index == 12
        homologous_radial1d_geometry.v_outer_boundary = (
                homologous_radial1d_geometry.v_outer[-2] - EPSILON_VELOCITY_SHIFT
        )
        assert homologous_radial1d_geometry.v_outer_boundary_index == 11

    def time_velocity_boundary(self):
        # testing the active cell boundaries when setting the boundaries

        homologous_radial1d_geometry = self.homologous_radial1d_geometry

        homologous_radial1d_geometry.v_inner_boundary = 7999 * u.km / u.s
        assert len(homologous_radial1d_geometry.v_inner_active) == len(
            homologous_radial1d_geometry.v_inner
        )

        homologous_radial1d_geometry.v_inner_boundary = 8001 * u.km / u.s
        assert len(homologous_radial1d_geometry.v_inner_active) == len(
            homologous_radial1d_geometry.v_inner
        )

        homologous_radial1d_geometry.v_inner_boundary = 9001 * u.km / u.s
        assert len(homologous_radial1d_geometry.v_inner_active) == (
                len(homologous_radial1d_geometry.v_inner) - 1
        )
        homologous_radial1d_geometry.v_inner_boundary = 9000 * u.km / u.s
        assert len(homologous_radial1d_geometry.v_inner_active) == (
                len(homologous_radial1d_geometry.v_inner) - 1
        )
