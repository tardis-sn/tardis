
from tardis.spectrum.base import TARDISSpectrum
from tardis.spectrum.formal_integral.source_function import SourceFunctionSolver
from tardis.spectrum.formal_integral.formal_integral_numba import NumbaFormalIntegrator
from tardis.spectrum.formal_integral.formal_integral_cuda import CudaFormalIntegrator


class FormalIntegralSolver:

    def __init__(self, formal_integral_configuration, spectrum_configuration, integrator_configuration):

        # will need configurations for
            # which method to use (numba, cuda, etc)
            # number of shells and points
            # spectrum configuration
        self.formal_integral_configuration = formal_integral_configuration
        self.spectrum_configuration = spectrum_configuration
        self.integrator_configuration = integrator_configuration

    def setup_optional(method):

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

        if method == 'numba':




        pass # TODO: Implement method to set up the integrator based on the configuration.

    def solve(self, sim_state, opacity_state, transport_state, atomic_data, levels):

        points = self.integrator_configuration.points
        interpolate_shells = self.integrator_configuration.interpolate_shells

        line_interaction_type = transport_state.line_interaction_type

        # TODO: determine configs for SourceFunctionSolver
        sourceFunction = SourceFunctionSolver(line_interaction_type=line_interaction_type)
        att_S_ul, Jredlu, Jbluelu, e_dot_u = sourceFunction.solve(sim_state, opacity_state, transport_state, atomic_data, levels)

        # interpolate if necessary
        # TODO: generalize the interpolation function
        # e.g. 
            # r_inner_i = interpolate(...)

        # TODO: compute the formal integral
        # 
        # L, I_nu_p = self.integrator.formal_integral()
        
        # TODO: make a TARDISSpectrum object from the results    
        frequency = 0  # Extract frequency from the results
        luminosity = 0  # Extract luminosity from the results
        spec = TARDISSpectrum(frequency, luminosity)
        return spec


    

