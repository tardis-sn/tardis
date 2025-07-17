import numpy as np
from scipy.interpolate import interp1d

from tardis.opacities.opacity_state import opacity_state_initialize
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry

from tardis.spectrum.base import TARDISSpectrum
from tardis.spectrum.formal_integral.base import check, interpolate_integrator_quantities
from tardis.spectrum.formal_integral.source_function import SourceFunctionSolver
from tardis.spectrum.formal_integral.formal_integral_numba import NumbaFormalIntegrator
from tardis.spectrum.formal_integral.formal_integral_cuda import CudaFormalIntegrator


class FormalIntegralSolver:

    def __init__(self, integrator_settings):

        # will need configurations for
            # which method to use (numba, cuda, etc)
            # number of shells and points
            # spectrum configuration
        # self.formal_integral_configuration = formal_integral_configuration
        # self.spectrum_configuration = spectrum_configuration
        self.integrator_settings = integrator_settings # TODO: add option to specify 'numba' or 'cuda'

    def setup(self, opacity_state, transport, plasma, macro_atom_state=None):

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
                transport.montecarlo_configuration
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

        return atomic_data, levels_index, opacity_state

    def setup_integrator(self, opacity_state, time_explosion, r_inner_i, r_outer_i):
        # TODO: move out - this has to happen after the interpolation, if any
        numba_radial_1d_geometry = NumbaRadial1DGeometry(
            r_inner_i,
            r_outer_i,
            r_inner_i
            / time_explosion.to("s").value,
            r_outer_i
            / time_explosion.to("s").value,
        )

        if self.integrator_settings.method == 'cuda':
            self.integrator = CudaFormalIntegrator(
                numba_radial_1d_geometry,
                time_explosion.cgs.value,
                opacity_state,
                self.integrator_settings.points,
            )
        else:
            self.integrator = NumbaFormalIntegrator(
                numba_radial_1d_geometry,
                time_explosion.cgs.value,
                opacity_state,
                self.integrator_settings.points
            )
    
    def solve(self, nu, simulation_state, opacity_state, transport, plasma, macro_atom_state=None):

        atomic_data, levels, opacity_state = self.setup(opacity_state, transport, plasma, macro_atom_state)
        transport_state = transport.transport_state

        points = self.integrator_settings.points
        interpolate_shells = self.integrator_settings.interpolate_shells
        line_interaction_type = transport.line_interaction_type

        sourceFunction = SourceFunctionSolver(line_interaction_type, atomic_data)
        res = sourceFunction.solve(simulation_state, opacity_state, transport_state, levels)

        att_S_ul, Jred_lu, Jblue_lu, e_dot_u = res.att_S_ul, res.Jred_lu, res.Jblue_lu, res.e_dot_u
        if interpolate_shells > 0: # TODO: fix up the interpolation
            (
                att_S_ul,
                Jred_lu,
                Jblue_lu,
                e_dot_u,
                r_inner_i,
                r_outer_i,
                tau_sobolevs_integ,
                electron_densities_integ
            ) = interpolate_integrator_quantities(
                att_S_ul, Jred_lu, Jblue_lu, e_dot_u,
                interpolate_shells,
                simulation_state, transport, opacity_state, plasma.electron_densities
            )
        else:
            r_inner_i = transport_state.geometry_state.r_inner
            r_outer_i = transport_state.geometry_state.r_outer
            tau_sobolevs_integ = opacity_state.tau_sobolev
            electron_densities_integ = opacity_state.electron_density

        att_S_ul = att_S_ul.flatten(order="F")
        Jred_lu = Jred_lu.flatten(order="F")
        Jblue_lu = Jblue_lu.flatten(order="F")

        self.setup_integrator(opacity_state, simulation_state.time_explosion, r_inner_i, r_outer_i)

        L, I_nu_p = self.integrator.formal_integral(
            simulation_state.t_inner,
            nu,
            nu.shape[0],
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            tau_sobolevs_integ,
            electron_densities_integ,
            points,
        )
        
        spec = TARDISSpectrum(nu, L)
        return spec
    

    # TODO: rewrite interpolate_integrator_quantities
    def interpolate_integrator_quantities(
        att_S_ul, Jredlu, Jbluelu, e_dot_u,
        interpolate_shells,
        simulation_state, transport, opacity_state, electron_densities
    ):
        """Interpolate the integrator quantities to interpolate_shells.

        Parameters
        ----------
        att_S_ul : np.ndarray
            attenuated source function for each line in each shell
        Jredlu : np.ndarray
            J estimator from the red end of the line from lower to upper level
        Jbluelu : np.ndarray
            J estimator from the blue end of the line from lower to upper level
        e_dot_u : np.ndarray
            Line estimator for the rate of energgity density absorption from lower to upper level
        interpolate_shells : int
            number of shells to interpolate to
        simulation_state : tardis.model.SimulationState
        transport : tardis.transport.montecarlo.MonteCarloTransportSolver
        opacity_state : OpacityStateNumba
        electron_densities : np.ndarray

        Returns
        -------
        tuple
            Interpolated values of att_S_ul, Jredlu, Jbluelu, and e_dot_u
        """

        mct_state = transport.transport_state
        
        nshells = interpolate_shells
        r_middle = (
            mct_state.geometry_state.r_inner + mct_state.geometry_state.r_outer
        ) / 2.0

        r_integ = np.linspace(
            mct_state.geometry_state.r_inner[0],
            mct_state.geometry_state.r_outer[-1],
            nshells,
        )
        r_inner_i = r_integ[:-1]
        r_outer_i = r_integ[1:]

        r_middle_integ = (r_integ[:-1] + r_integ[1:]) / 2.0

        electron_densities_integ = interp1d(
            r_middle,
            electron_densities.iloc[
                simulation_state.geometry.v_inner_boundary_index : simulation_state.geometry.v_outer_boundary_index
            ],
            fill_value="extrapolate",
            kind="nearest",
        )(r_middle_integ)
        # Assume tau_sobolevs to be constant within a shell
        # (as in the MC simulation)
        tau_sobolevs_integ = interp1d(
            r_middle,
            opacity_state.tau_sobolev[
                :,
                simulation_state.geometry.v_inner_boundary_index : simulation_state.geometry.v_outer_boundary_index,
            ],
            fill_value="extrapolate",
            kind="nearest",
        )(r_middle_integ)
        att_S_ul = interp1d(r_middle, att_S_ul, fill_value="extrapolate")(
            r_middle_integ
        )
        Jredlu = interp1d(r_middle, Jredlu, fill_value="extrapolate")(
            r_middle_integ
        )
        Jbluelu = interp1d(r_middle, Jbluelu, fill_value="extrapolate")(
            r_middle_integ
        )
        e_dot_u = interp1d(r_middle, e_dot_u, fill_value="extrapolate")(
            r_middle_integ
        )

        # Set negative values from the extrapolation to zero
        att_S_ul = att_S_ul.clip(0.0)
        Jbluelu = Jbluelu.clip(0.0)
        Jredlu = Jredlu.clip(0.0)
        e_dot_u = e_dot_u.clip(0.0)
        return att_S_ul, Jredlu, Jbluelu, e_dot_u, r_inner_i, r_outer_i, tau_sobolevs_integ, electron_densities_integ




    

