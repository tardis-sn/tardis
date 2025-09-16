import logging

import numpy as np
from astropy import units as u
from scipy.interpolate import interp1d

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.spectrum.base import TARDISSpectrum
from tardis.spectrum.formal_integral.base import check_formal_integral_requirements
from tardis.spectrum.formal_integral.formal_integral_cuda import (
    CudaFormalIntegrator,
)
from tardis.spectrum.formal_integral.formal_integral_numba import (
    NumbaFormalIntegrator,
)
from tardis.spectrum.formal_integral.source_function import SourceFunctionSolver

logger = logging.getLogger(__name__)


class FormalIntegralSolver:
    """
    Formal integral solver for TARDIS spectra.

    Attributes
    ----------
    points : int
        Number of points for the formal integral calculation
    interpolate_shells : int
        Number of shells to interpolate to. If > 0, interpolation is performed
        to the specified number of shells. If <= 0, no interpolation is performed
        and original shell structure is used.
    method : str or None
        Method to use for the formal integral solver ('numba' or 'cuda')
    """

    points: int
    interpolate_shells: int
    method: str | None

    def __init__(
        self, points: int, interpolate_shells: int, method: str | None = None
    ) -> None:
        """
        Initialize the formal integral solver.

        Parameters
        ----------
        points : int
            Number of points for the formal integral calculation
        interpolate_shells : int
            Number of shells to interpolate to. If > 0, interpolation is performed
            to the specified number of shells. If <= 0, no interpolation is performed
            and original shell structure is used.
        method : str, optional
            Method to use for the formal integral solver ('numba' or 'cuda').
            If None, will be determined based on GPU availability.
        """
        self.points = points
        self.interpolate_shells = interpolate_shells
        self.method = method

    def setup(
        self,
        transport,
        opacity_state=None,
        macro_atom_state=None,
    ) -> object:
        """
        Prepare the necessary data for the formal integral solver.

        Parameters
        ----------
        transport : tardis.transport.montecarlo.MontecarloTransport
            The transport configuration object
        opacity_state : tardis.opacities.opacity_state.OpacityState
            The regular (non-numba) opacity state object to be converted
        macro_atom_state : tardis.opacities.macro_atom.macroatom_state.MacroAtomState, optional
            The macro atom state object

        Returns
        -------
        opacity_state_numba : tardis.opacities.opacity_state.OpacityStateNumba
            The opacity state converted to numba format
        """
        if opacity_state is None or macro_atom_state is None:
            raise NotImplementedError(
                "This functionality does not work anymore. Both opacity_state and macro_atom_state must be provided."
            )

        self.montecarlo_configuration = transport.montecarlo_configuration

        if self.method in [
            None,
            "numba",
            "cuda",
        ]:  # TODO: better way to handle this
            # use GPU if available
            if transport.use_gpu:
                self.method = "cuda"
            else:
                self.method = "numba"
        else:
            logger.warning(
                "Computing formal integral via the %s method isn't supported. "
                "Please run with config option numba or cuda. "
                "Defaulting to numba implementation",
                self.method,
            )
            self.method = "numba"

        if opacity_state and macro_atom_state:
            opacity_state_numba = opacity_state.to_numba(
                macro_atom_state,
                transport.line_interaction_type,
            )

        return opacity_state_numba

    def setup_integrator(
        self,
        opacity_state_numba,
        time_explosion: u.Quantity,
        r_inner: np.ndarray,
        r_outer: np.ndarray,
    ) -> None:
        """
        Set up the integrator depending on the choice of method.

        Parameters
        ----------
        opacity_state_numba : tardis.opacities.opacity_state.OpacityStateNumba
            The opacity state (numba format) to use for the formal integral
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
                opacity_state_numba,
                self.points,
            )
        else:
            self.integrator = NumbaFormalIntegrator(
                numba_radial_1d_geometry,
                time_explosion.cgs.value,
                opacity_state_numba,
                self.points,
            )

    def solve(
        self,
        frequencies: u.Quantity,
        simulation_state,
        transport_solver,
        opacity_state,
        atomic_data,
        electron_densities,
        macro_atom_state=None,
    ) -> TARDISSpectrum:
        """
        Solve the formal integral.

        Parameters
        ----------
        frequencies : u.Quantity
            The frequency grid for the formal integral
        simulation_state : tardis.model.SimulationState
            State which holds information about each shell
        transport_solver : tardis.transport.montecarlo.MonteCarloTransportSolver
            The transport solver
        opacity_state : tardis.opacities.opacity_state.OpacityState
            Regular (non-numba) opacity state; will be converted to numba via `setup`
        atomic_data : tardis.atomic.AtomicData
            Atomic data containing atomic properties
        electron_densities : pd.Series
            Electron densities for each shell
        macro_atom_state : tardis.opacities.macro_atom.macroatom_state.MacroAtomState, optional
            State of the macro atom (required for converting opacity_state to numba)

        Returns
        -------
        TARDISSpectrum
            The formal integral spectrum
        """
        # check objects and configs
        check_formal_integral_requirements(simulation_state, opacity_state, transport_solver)

        # Convert to numba opacity state for source function and integrator
        opacity_state_numba = self.setup(
            transport_solver, opacity_state, macro_atom_state
        )
        transport_state = transport_solver.transport_state

        points = self.points
        interpolate_shells = self.interpolate_shells
        line_interaction_type = transport_solver.line_interaction_type

        if interpolate_shells == 0:  # Default Value
            interpolate_shells = max(2 * simulation_state.no_of_shells, 80)
            logger.warning(
                "The number of interpolate_shells was not "
                "specified. The value was set to %d.",
                interpolate_shells
            )
        self.interpolate_shells = interpolate_shells

        source_function_solver = SourceFunctionSolver(line_interaction_type)
        source_function_state = source_function_solver.solve(
            simulation_state, opacity_state_numba, transport_state, atomic_data
        )

        # Generate interpolated radii if needed
        mct_state = transport_solver.transport_state
        if interpolate_shells > 0:
            radius_interpolated = np.linspace(
                mct_state.geometry_state.r_inner[0],
                mct_state.geometry_state.r_outer[-1],
                interpolate_shells,
            )
            r_inner_interpolated = radius_interpolated[:-1]
            r_outer_interpolated = radius_interpolated[1:]
        elif interpolate_shells <= 0:
            # Use original radii values when interpolate_shells < 0
            r_inner_interpolated = mct_state.geometry_state.r_inner
            r_outer_interpolated = mct_state.geometry_state.r_outer

        (
            att_S_ul_interpolated,
            Jred_lu_interpolated,
            Jblue_lu_interpolated,
            r_inner_interpolated,
            r_outer_interpolated,
            tau_sobolevs_interpolated,
            electron_densities_interpolated,
        ) = self.interpolate_integrator_quantities(
            mct_state.geometry_state.r_inner,
            mct_state.geometry_state.r_outer,
            r_inner_interpolated,
            r_outer_interpolated,
            source_function_state,
            simulation_state,
            opacity_state,
            electron_densities,
        )

        att_S_ul_interpolated = att_S_ul_interpolated.flatten(order="F")
        Jred_lu_interpolated = Jred_lu_interpolated.flatten(order="F")
        Jblue_lu_interpolated = Jblue_lu_interpolated.flatten(order="F")

        self.setup_integrator(
            opacity_state_numba,
            simulation_state.time_explosion,
            r_inner_interpolated,
            r_outer_interpolated,
        )

        luminosity_densities, intensities_nu_p = self.integrator.formal_integral(
            simulation_state.t_inner,
            frequencies,
            att_S_ul_interpolated,
            Jred_lu_interpolated,
            Jblue_lu_interpolated,
            tau_sobolevs_interpolated,
            electron_densities_interpolated,
            points,
        )

        luminosity_densities = np.array(luminosity_densities, dtype=np.float64)
        delta_frequency = frequencies[1] - frequencies[0]

        assert np.allclose(
            frequencies.diff(), delta_frequency, atol=0, rtol=1e-12
        ), "Frequency grid must be uniform"

        luminosity = u.Quantity(luminosity_densities, "erg/s/Hz") * delta_frequency

        frequencies = frequencies.to("Hz", u.spectral())

        # Ugly hack to convert to 'bin edges'
        frequencies = u.Quantity(
            np.concatenate(
                [
                    frequencies.value,
                    [frequencies.value[-1] + np.diff(frequencies.value)[-1]],
                ]
            ),
            frequencies.unit,
        )

        return TARDISSpectrum(frequencies, luminosity)

    def interpolate_integrator_quantities(
        self,
        r_inner_original: np.ndarray,
        r_outer_original: np.ndarray,
        r_inner_interpolated: np.ndarray,
        r_outer_interpolated: np.ndarray,
        source_function_state,
        simulation_state,
        opacity_state,
        electron_densities,
    ) -> tuple[
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
    ]:
        """
        Interpolate the integrator quantities to interpolate_shells.

        Parameters
        ----------
        r_inner_original : np.ndarray
            Original inner radii of the shells
        r_outer_original : np.ndarray
            Original outer radii of the shells
        r_inner_interpolated : np.ndarray
            Pre-computed inner radii for interpolation
        r_outer_interpolated : np.ndarray
            Pre-computed outer radii for interpolation
        source_function_state : tardis.spectrum.formal_integral.source_function.SourceFunctionState
            Data class that holds the computed source function values which will be interpolated
        simulation_state : tardis.model.SimulationState
            The simulation state object
        opacity_state : tardis.opacities.opacity_state.OpacityState
            The opacity state object (regular, non-numba)
        electron_densities : pd.Series
            Electron densities for each shell

        Returns
        -------
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
            Interpolated values of att_S_ul_interpolated, Jred_lu_interpolated, Jblue_lu_interpolated, r_inner_interpolated, r_outer_interpolated, tau_sobolevs_interpolated, and electron_densities_interpolated
        """
        # Extract values from source function state
        att_S_ul = source_function_state.att_S_ul
        Jred_lu = source_function_state.Jred_lu
        Jblue_lu = source_function_state.Jblue_lu

        r_middle_original = (r_inner_original + r_outer_original) / 2.0

        r_middle_interpolated = (r_inner_interpolated + r_outer_interpolated) / 2.0

        electron_densities_interpolated = interp1d(
            r_middle_original,
            electron_densities.iloc[
                simulation_state.geometry.v_inner_boundary_index : simulation_state.geometry.v_outer_boundary_index
            ],
            fill_value="extrapolate",  # type: ignore[arg-type]
            kind="nearest",
        )(r_middle_interpolated)
        # Assume tau_sobolevs to be constant within a shell
        # (as in the MC simulation)
        v_inner_boundary_index = simulation_state.geometry.v_inner_boundary_index
        v_outer_boundary_index = simulation_state.geometry.v_outer_boundary_index
        tau_sobolevs_interpolated = interp1d(
            r_middle_original,
            opacity_state.tau_sobolev.values[
                :, v_inner_boundary_index:v_outer_boundary_index
            ],
            fill_value="extrapolate",  # type: ignore[arg-type]
            kind="nearest",
        )(r_middle_interpolated)
        att_S_ul_interpolated = interp1d(
            r_middle_original,
            att_S_ul,
            fill_value="extrapolate",  # type: ignore[arg-type]
        )(r_middle_interpolated)
        Jred_lu_interpolated = interp1d(
            r_middle_original,
            Jred_lu,
            fill_value="extrapolate",  # type: ignore[arg-type]
        )(r_middle_interpolated)
        Jblue_lu_interpolated = interp1d(
            r_middle_original,
            Jblue_lu,
            fill_value="extrapolate",  # type: ignore[arg-type]
        )(r_middle_interpolated)

        # Set negative values from the extrapolation to zero
        att_S_ul_interpolated = att_S_ul_interpolated.clip(0.0)
        Jblue_lu_interpolated = Jblue_lu_interpolated.clip(0.0)
        Jred_lu_interpolated = Jred_lu_interpolated.clip(0.0)
        return (
            att_S_ul_interpolated,
            Jred_lu_interpolated,
            Jblue_lu_interpolated,
            r_inner_interpolated,
            r_outer_interpolated,
            tau_sobolevs_interpolated,
            electron_densities_interpolated,
        )
