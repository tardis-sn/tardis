import warnings

import numpy as np
from astropy import units as u
from scipy.interpolate import interp1d

from tardis.opacities.opacity_state import opacity_state_initialize
from tardis.spectrum.formal_integral.formal_integral_cuda import (
    CudaFormalIntegrator,
)
from tardis.spectrum.formal_integral.formal_integral_numba import (
    NumbaFormalIntegrator,
    calculate_p_values,
    trapezoid_integration,
)
from tardis.spectrum.formal_integral.base import (
    check,
    interpolate_integrator_quantities,
)
from tardis.spectrum.formal_integral.source_function import SourceFunctionSolver
from tardis.spectrum.spectrum import TARDISSpectrum


class FormalIntegrator:
    """
    Class containing the formal integrator.

    If there is a NVIDIA CUDA GPU available,
    the formal integral will automatically run
    on it. If multiple GPUs are available, it will
    choose the first one that it sees. You can
    read more about selecting different GPUs on
    Numba's CUDA documentation.

    Parameters
    ----------
    model : tardis.model.SimulationState
    plasma : tardis.plasma.BasePlasma
    transport : tardis.transport.montecarlo.MontecarloTransport
    points : int64
    """

    def __init__(
        self,
        simulation_state,
        plasma,
        transport,
        opacity_state=None,
        macro_atom_state=None,
        points=1000,
    ):
        self.simulation_state = simulation_state
        self.transport = transport
        self.points = points
        if transport:
            self.montecarlo_configuration = (
                self.transport.montecarlo_configuration
            )
        if plasma and opacity_state and macro_atom_state:
            self.opacity_state = opacity_state.to_numba(
                macro_atom_state,
                transport.line_interaction_type,
            )
            self.atomic_data = plasma.atomic_data
            self.plasma = plasma
            self.levels_index = plasma.levels
        elif plasma:
            self.opacity_state = opacity_state_initialize(
                plasma,
                transport.line_interaction_type,
                self.montecarlo_configuration.DISABLE_LINE_SCATTERING,
            )
            self.atomic_data = plasma.atomic_data
            self.plasma = plasma
            self.levels_index = plasma.levels
        else:
            self.opacity_state = None

    def generate_numba_objects(self):
        """
        Instantiate the numba interface objects
        needed for computing the formal integral
        """
        from tardis.model.geometry.radial1d import NumbaRadial1DGeometry

        self.numba_radial_1d_geometry = NumbaRadial1DGeometry(
            self.transport.r_inner_i,
            self.transport.r_outer_i,
            self.transport.r_inner_i
            / self.simulation_state.time_explosion.to("s").value,
            self.transport.r_outer_i
            / self.simulation_state.time_explosion.to("s").value,
        )
        if self.opacity_state is None:
            self.opacity_state = opacity_state_initialize(
                self.plasma,
                self.transport.line_interaction_type,
                self.montecarlo_configuration.DISABLE_LINE_SCATTERING,
            )
        if self.transport.use_gpu:
            self.integrator = CudaFormalIntegrator(
                self.numba_radial_1d_geometry,
                self.simulation_state.time_explosion.cgs.value,
                self.opacity_state,
                self.points,
            )
        else:
            self.integrator = NumbaFormalIntegrator(
                self.numba_radial_1d_geometry,
                self.simulation_state.time_explosion.cgs.value,
                self.opacity_state,
                self.points,
            )

    def calculate_spectrum(
        self, frequency, points=None, interpolate_shells=0, raises=True
    ):
        # Very crude implementation
        # The c extension needs bin centers (or something similar)
        # while TARDISSpectrum needs bin edges
        check(self.simulation_state, self.plasma, self.transport, raises=raises)
        N = points or self.points
        if interpolate_shells == 0:  # Default Value
            interpolate_shells = max(2 * self.simulation_state.no_of_shells, 80)
            warnings.warn(
                "The number of interpolate_shells was not "
                f"specified. The value was set to {interpolate_shells}."
            )
        self.interpolate_shells = interpolate_shells
        frequency = frequency.to("Hz", u.spectral())

        luminosity = u.Quantity(self.formal_integral(frequency, N), "erg") * (
            frequency[1] - frequency[0]
        )

        # Ugly hack to convert to 'bin edges'
        frequency = u.Quantity(
            np.concatenate(
                [
                    frequency.value,
                    [frequency.value[-1] + np.diff(frequency.value)[-1]],
                ]
            ),
            frequency.unit,
        )

        return TARDISSpectrum(frequency, luminosity)

    def formal_integral(self, nu, N):
        """
        Do the formal integral with the numba routines
        """
        # TODO: get rid of storage later on

        transport_state = self.transport.transport_state

        source_function_solver = SourceFunctionSolver(
            self.transport.line_interaction_type, self.plasma.atomic_data
        )
        source_function_state = source_function_solver.solve(
            self.simulation_state,
            self.opacity_state,
            transport_state,
            self.plasma.levels,
        )

        if self.interpolate_shells > 0:
            (
                att_S_ul,
                Jred_lu,
                Jblue_lu,
                e_dot_u,
            ) = interpolate_integrator_quantities(
                source_function_state,
                self.interpolate_shells,
                self.simulation_state,
                self.transport,
                self.opacity_state,
                self.plasma.electron_densities,
            )
        else:
            self.transport.r_inner_i = transport_state.geometry_state.r_inner
            self.transport.r_outer_i = transport_state.geometry_state.r_outer
            self.transport.tau_sobolevs_integ = self.opacity_state.tau_sobolev
            self.transport.electron_densities_integ = (
                self.opacity_state.electron_density
            )

        att_S_ul = att_S_ul.flatten(order="F")
        Jred_lu = Jred_lu.flatten(order="F")
        Jblue_lu = Jblue_lu.flatten(order="F")

        self.generate_numba_objects()
        L, I_nu_p = self.integrator.formal_integral(
            self.simulation_state.t_inner,
            nu,
            nu.shape[0],
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            self.transport.tau_sobolevs_integ,
            self.transport.electron_densities_integ,
            N,
        )
        R_max = self.transport.r_outer_i[-1]
        ps = calculate_p_values(R_max, N)[None, :]
        I_nu_p[:, 1:] /= ps[:, 1:]
        self.transport.I_nu_p = I_nu_p
        self.transport.p_rays = ps

        I_nu = self.transport.I_nu_p * ps
        L_test = np.array(
            [
                8
                * np.pi
                * np.pi
                * trapezoid_integration((I_nu)[i, :], R_max / N)
                for i in range(nu.shape[0])
            ]
        )
        error = np.max(np.abs((L_test - L) / L))
        assert error < 1e-7, (
            f"Incorrect I_nu_p values, max relative difference:{error}"
        )

        return np.array(L, np.float64)
