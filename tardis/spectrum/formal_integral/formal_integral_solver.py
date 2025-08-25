import numpy as np
from astropy import units as u
import warnings
from scipy.interpolate import interp1d
import logging

from tardis.opacities.opacity_state import opacity_state_initialize
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry

from tardis.spectrum.base import TARDISSpectrum
from tardis.spectrum.formal_integral.base import check
from tardis.spectrum.formal_integral.source_function import SourceFunctionSolver
from tardis.spectrum.formal_integral.formal_integral_numba import (
    NumbaFormalIntegrator,
)
from tardis.spectrum.formal_integral.formal_integral_cuda import (
    CudaFormalIntegrator,
)

logger = logging.getLogger(__name__)


class FormalIntegralSolver:
    def __init__(self, integrator_settings):
        """
        Initialize the formal integral solver.

        Parameters
        ----------
        integrator_settings : IntegratorSettings
            The settings to use for the integrator, such as:
                points (int): Number of points
                interpolate_shells (int): Number of shells to interpolate to
                method (str): Method to use for the formal integral solver ('numba' or 'cuda')
        """
        self.integrator_settings = integrator_settings

        # check if the method was set in the configuration
        try:
            self.method = self.integrator_settings.integrated.method
        except AttributeError:
            self.method = None

    def setup(
        self, transport, plasma, opacity_state=None, macro_atom_state=None
    ):
        """
        Prepares the necessary data for the formal integral solver.

        Parameters
        ----------
        plasma : tardis.plasma.BasePlasma
        transport : tardis.transport.montecarlo.MontecarloTransport

        Returns
        -------
        atomic_data : AtomicData
        levels_index : np.ndarray
        opacity_state : OpacityStateNumba
        """

        self.montecarlo_configuration = transport.montecarlo_configuration

        # check method selection
        if self.method in [ # TODO: better way to handle this
            "numba",
            "cuda",
        ]:
            pass
        elif self.method is None:
            logger.warning(
                f"The formal integral implementation was not specified. "
                "Please run with config option numba or cuda"
                "Defaulting to numba implementation"
            )
            self.method = "numba"
        else:
            logger.warning(
                f"Computing formal integral via the {self.method} method isn't supported"
                "Please run with config option numba or cuda"
                "Defaulting to numba implementation"
            )
            self.method = "numba"

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

        return atomic_data, opacity_state

    def setup_integrator(self, opacity_state, time_explosion, r_inner, r_outer):
        """
        Setup the integrator depending on the choice of method

        Parameters
        ----------
        opacity_state : tardis.opacities.opacity_state.OpacityStateNumba
            The opacity state to use for the formal integral
        time_explosion : u.Quantity
            The time of the explosion
        r_inner : np.ndarray
            The inner radii of the shells
        r_outer : np.ndarray
            The outer radii of the shells
        """

        numba_radial_1d_geometry = NumbaRadial1DGeometry(
            r_inner,
            r_outer,
            r_inner / time_explosion.to("s").value,
            r_outer / time_explosion.to("s").value,
        )

        if self.method == "cuda":
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
                self.integrator_settings.points,
            )

    def solve(
        self,
        nu,
        simulation_state,
        transport,
        plasma,
        opacity_state=None,
        macro_atom_state=None,
    ):
        """
        Solve the formal integral

        Parameters
        ----------
        nu : u.Quantity
            The frequency grid for the formal integral.
        simulation_state : tardis.model.SimulationState
            State which hold information about each shell
        transport : tardis.transport.montecarlo.MonteCarloTransportSolver
        plasma : tardis.plasma.BasePlasma
        opacity_state : tardis.opacities.opacity_state.OpacityState or None
            State of the line opacities
        macro_atom_state : tardis.opacities.macro_atom.macroatom_state.MacroAtomState or None
            State of the macro atom

        Returns
        -------
        TARDISSpectrum
            the formal integral spectrum
        """

        atomic_data, opacity_state = self.setup(
            transport, plasma, opacity_state, macro_atom_state
        )
        transport_state = transport.transport_state

        points = self.integrator_settings.points
        interpolate_shells = self.integrator_settings.interpolate_shells

        if interpolate_shells == 0:  # Default Value
            interpolate_shells = max(2 * simulation_state.no_of_shells, 80)
            warnings.warn(
                "The number of interpolate_shells was not "
                f"specified. The value was set to {interpolate_shells}."
            )
        self.integrator_settings.interpolate_shells = interpolate_shells

        line_interaction_type = transport.line_interaction_type

        source_function_solver = SourceFunctionSolver(line_interaction_type)
        source_function_state = source_function_solver.solve(
            simulation_state, opacity_state, transport_state, atomic_data
        )

        (
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            _,  # e_dot_u is not used
            r_inner_itp,
            r_outer_itp,
            tau_sobolevs_integ,
            electron_densities_integ,
        ) = self.get_interpolated_quantities(
            source_function_state,
            interpolate_shells,
            simulation_state,
            transport,
            opacity_state,
            plasma,
        )
        att_S_ul = att_S_ul.flatten(order="F")
        Jred_lu = Jred_lu.flatten(order="F")
        Jblue_lu = Jblue_lu.flatten(order="F")

        self.setup_integrator(
            opacity_state,
            simulation_state.time_explosion,
            r_inner_itp,
            r_outer_itp,
        )

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

        L = np.array(L, dtype=np.float64)
        luminosity = u.Quantity(L, "erg") * (nu[1] - nu[0])

        frequency = nu.to("Hz", u.spectral())

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

    # TODO: rewrite interpolate_integrator_quantities
    def interpolate_integrator_quantities(
        self,
        att_S_ul,
        Jredlu,
        Jbluelu,
        e_dot_u,
        interpolate_shells,
        simulation_state,
        transport,
        opacity_state,
        electron_densities,
    ):
        """
        Interpolate the integrator quantities to interpolate_shells.

        Parameters
        ----------
        source_function_state : tardis.spectrum.formal_integral.source_function.SourceFunctionState
            Data class that hold the computed source function values
        interpolate_shells : int
            number of shells to interpolate to
        simulation_state : tardis.model.SimulationState
        transport : tardis.transport.montecarlo.MonteCarloTransportSolver
        opacity_state : tardis.opacities.opacity_state.OpacityStateNumba
        electron_densities : np.ndarray

        Returns
        -------
        tuple
            Interpolated values of att_S_ul, Jredlu, Jbluelu, e_dot_u, r_inner_i, r_outer_i, tau_sobolevs_integ, and electron_densities_integ
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
        return (
            att_S_ul,
            Jredlu,
            Jbluelu,
            e_dot_u,
            r_inner_i,
            r_outer_i,
            tau_sobolevs_integ,
            electron_densities_integ,
        )

    def get_interpolated_quantities(
        self,
        source_function_state,
        interpolate_shells,
        simulation_state,
        transport,
        opacity_state,
        plasma,
    ):
        """
        If needed, interpolate the quantities from the source function state, and prepare the results for use in the formal integral.

        Parameters
        ----------
        source_function_state : tardis.spectrum.formal_integral.source_function.SourceFunctionState
            Data class that hold the computed source function values which will be interpolated, if needed
        interpolate_shells : int
            The number of shells to interpolate to.
        simulation_state : tardis.model.SimulationState
        transport : tardis.transport.montecarlo.MonteCarloTransportSolver
        opacity_state : tardis.opacities.opacity_state.OpacityStateNumba
        plasma : tardis.plasma.Plasma

        Returns
        -------
        tuple
            (possibly interpolated) att_S_ul, Jred_lu, Jblue_lu, e_dot_u, r_inner, r_outer, tau_sobolevs, electron_densities
        """

        att_S_ul, Jred_lu, Jblue_lu, e_dot_u = (
            source_function_state.att_S_ul,
            source_function_state.Jred_lu,
            source_function_state.Jblue_lu,
            source_function_state.e_dot_u,
        )

        # interpolate, if not use existing values
        if interpolate_shells > 0:
            (
                att_S_ul,
                Jred_lu,
                Jblue_lu,
                e_dot_u,
                r_inner_i,
                r_outer_i,
                tau_sobolevs_integ,
                electron_densities_integ,
            ) = self.interpolate_integrator_quantities(
                att_S_ul,
                Jred_lu,
                Jblue_lu,
                e_dot_u,
                interpolate_shells,
                simulation_state,
                transport,
                opacity_state,
                plasma.electron_densities,
            )
        else:
            r_inner_i = transport.transport_state.geometry_state.r_inner
            r_outer_i = transport.transport_state.geometry_state.r_outer
            tau_sobolevs_integ = opacity_state.tau_sobolev
            electron_densities_integ = opacity_state.electron_density

        return (
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            e_dot_u,
            r_inner_i,
            r_outer_i,
            tau_sobolevs_integ,
            electron_densities_integ,
        )
