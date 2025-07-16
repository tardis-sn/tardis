from tardis.opacities.opacity_state import opacity_state_initialize
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry

from tardis.spectrum.base import TARDISSpectrum
from tardis.spectrum.formal_integral.base import check, interpolate_integrator_quantities
from tardis.spectrum.formal_integral.source_function import SourceFunctionSolver
from tardis.spectrum.formal_integral.formal_integral_numba import NumbaFormalIntegrator
from tardis.spectrum.formal_integral.formal_integral_cuda import CudaFormalIntegrator


class FormalIntegralSolver:

    def __init__(self, integrator_configuration):

        # will need configurations for
            # which method to use (numba, cuda, etc)
            # number of shells and points
            # spectrum configuration
        # self.formal_integral_configuration = formal_integral_configuration
        # self.spectrum_configuration = spectrum_configuration
        self.integrator_configuration = integrator_configuration # TODO: add option to specify 'numba' or 'cuda'
        # self.montecarlo_configuration = montecarlo_configuration

    def setup(self, simulation_state, opacity_state, transport, plasma, macro_atom_state=None):

        """
        Set up the integrator depending on the method specified in the configuration.

        Parameters
        ----------
        method : str
            The method to use for the formal integral solver, e.g., 'numba', 'cuda'.

        Returns
        -------
        integrator : FormalIntegrator
            An instance of the appropriate formal integrator.
        """
        # atomic_data = plasma.atomic_data
        # levels = plasma.levels

        if transport:
            self.montecarlo_configuration = (
                self.transport.montecarlo_configuration
            )

        # TODO warn if no plasma
        if opacity_state and macro_atom_state:
            opacity_state = opacity_state.to_numba(
                macro_atom_state,
                transport.line_interaction_type,
            )
        else:
            opacity_state = opacity_state_initialize(
                plasma,
                transport.line_interaction_type,
                self.montecarlo_configuration.DISABLE_LINE_SCATTERING,
            )
        atomic_data = plasma.atomic_data
        levels_index = plasma.levels

        numba_radial_1d_geometry = NumbaRadial1DGeometry(
            transport.r_inner_i,
            transport.r_outer_i,
            transport.r_inner_i
            / simulation_state.time_explosion.to("s").value,
            transport.r_outer_i
            / simulation_state.time_explosion.to("s").value,
        )

        if self.integrator_configuration.method == 'cuda':
            self.integrator = CudaFormalIntegrator(
                numba_radial_1d_geometry,
                simulation_state.time_explosion.cgs.value,
                opacity_state,
                self.integrator_configuration.points,
            )
        else:
            self.integrator = NumbaFormalIntegrator(
                numba_radial_1d_geometry,
                simulation_state.time_explosion.cgs.value,
                opacity_state,
                self.integrator_configuration.points
            )

        return atomic_data, levels_index, opacity_state

    def solve(self, nu, simulation_state, opacity_state, transport, plasma):

        atomic_data, levels, opacity_state = self.setup(simulation_state, opacity_state, transport, plasma)
        transport_state = transport.transport_state

        points = self.integrator_configuration.points
        interpolate_shells = self.integrator_configuration.interpolate_shells
        line_interaction_type = transport_state.line_interaction_type

        sourceFunction = SourceFunctionSolver(line_interaction_type=line_interaction_type)
        res = sourceFunction.solve(simulation_state, opacity_state, transport_state, atomic_data, levels)

        att_S_ul, Jred_lu, Jblue_lu, e_dot_u = res.att_S_ul, res.Jred_lu, res.Jblue_lu, res.e_dot_u
        if self.interpolate_shells > 0: # TODO: fix up the interpolation
            (
                att_S_ul,
                Jred_lu,
                Jblue_lu,
                e_dot_u,
            ) = interpolate_integrator_quantities(
                att_S_ul, Jred_lu, Jblue_lu, e_dot_u,
                interpolate_shells,
                simulation_state, transport, opacity_state, plasma.electron_densities
            )
        else:
            # TODO: check jax version where I sorted all the bits
            self.transport.r_inner_i = ( #TODO: fix so this doesnt happen!!! no assigning stuff to the transport
                transport_state.geometry_state.r_inner
            )
            self.transport.r_outer_i = (
                transport_state.geometry_state.r_outer
            )
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
            points,
        )
        
        spec = TARDISSpectrum(nu, L)
        return spec
    

    # TODO: rewrite interpolate_integrator_quantities
    def interpolate():




    

