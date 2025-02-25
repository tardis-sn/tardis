import logging

import numpy as np
import pandas as pd
from astropy import units as u
from scipy.interpolate import interp1d

from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.simulation.convergence import ConvergenceSolver
from tardis.workflows.simple_tardis_workflow import SimpleTARDISWorkflow
from tardis.workflows.util import get_tau_integ

# logging support
logger = logging.getLogger(__name__)

# TODO:
# Option to estimate initial v_inner from electron opacity
# Add option for desired optical depth
# Think about csvy vs other formats for specifying grid
# Handle non-explicit formats when going out of the simulation


class InnerVelocitySolverWorkflow(SimpleTARDISWorkflow):
    TAU_TARGET = np.log(2.0 / 3.0)

    def __init__(self, configuration, mean_optical_depth="rosseland", tau=None):
        super().__init__(configuration)
        self.mean_optical_depth = mean_optical_depth.lower()

        self.convergence_solvers["v_inner_boundary"] = ConvergenceSolver(
            self.convergence_strategy.v_inner_boundary
        )

        # Need to compute the opacity state on init to get the optical depths
        # for the first inner boundary calculation.
        self.opacity_states = self.solve_opacity()

        if tau is not None:
            self.TAU_TARGET = np.log(tau)

        self.iterations_w = np.full((self.total_iterations, self.simulation_state.no_of_shells), np.nan)
        self.iterations_t_rad = np.full((self.total_iterations, self.simulation_state.no_of_shells), np.nan) * u.K
        self.iterations_electron_densities = np.full((self.total_iterations, self.simulation_state.no_of_shells), np.nan)
        self.iterations_t_inner = np.full(self.total_iterations, np.nan) * u.K
        self.iterations_v_inner_boundary = np.full(self.total_iterations, np.nan) * u.cm / u.s
        self.iterations_mean_optical_depth = np.full((self.total_iterations, self.simulation_state.no_of_shells), np.nan)

        initial_v_inner = self.estimate_v_inner()

        self.simulation_state.geometry.v_inner_boundary = initial_v_inner
        self.simulation_state.blackbody_packet_source.radius = (
            self.simulation_state.r_inner[0]
        )

    def store_plasma_state(self, i, t_radiative, dilution_factor, electron_densities, t_inner, v_inner_boundary, tau_integ):
        """Store current plasma information, 
        including the velocity of the inner boundary 
        and the Rosseland mean optical depth,
        used in iterated i.

        Parameters
        ----------
        i : int
            current iteration index (0 for the first)
        t_rad : astropy.units.Quantity
            radiation temperature
        dilution_factor : np.ndarray
            dilution factor
        electron_densities : np.ndarray
            electron density
        t_inner : astropy.units.Quantity
            temperature of inner boundary
        v_inner_boundary : astropy.units.Quantity
            velocity of inner boundary
        tau_integ : np.ndarray
            Rosseland mean optical depth
        """
        self.iterations_t_rad[i, -len(t_radiative):] = t_radiative
        self.iterations_w[i, -len(dilution_factor):] = dilution_factor
        self.iterations_electron_densities[i, -len(electron_densities):] = electron_densities.values
        self.iterations_t_inner[i] = t_inner
        self.iterations_v_inner_boundary[i] = v_inner_boundary
        self.iterations_mean_optical_depth[i,-len(tau_integ):] = tau_integ

    def reshape_store_plasma_state(self, executed_iterations):
        """Reshapes the storage arrays in case convergence was reached before
        all specified iterations were executed.

        Parameters
        ----------
        executed_iterations : int
            iteration index, i.e. number of iterations executed minus one!
        """
        self.iterations_t_rad = self.iterations_t_rad[
            : executed_iterations + 1, :
        ]
        self.iterations_w = self.iterations_w[: executed_iterations + 1, :]
        self.iterations_electron_densities = self.iterations_electron_densities[
            : executed_iterations + 1, :
        ]
        self.iterations_t_inner = self.iterations_t_inner[
            : executed_iterations + 1
        ]
        self.iterations_v_inner_boundary = self.iterations_v_inner_boundary[
            : executed_iterations + 1
        ]
        self.iterations_mean_optical_depth = self.iterations_mean_optical_depth[
            : executed_iterations + 1, :
        ]

    def estimate_v_inner(self):
        """
        Compute the Rosseland Mean Optical Depth,
        Estimate location where v_inner makes t=2/3 (or target)
        Extrapolate with exponential fits

        Need some way to return and inspect the optical depths for later logging
        """

        self.tau_integ = np.log(
            get_tau_integ(
                self.plasma_solver,
                self.opacity_states["opacity_state"],
                self.simulation_state,
            )[self.mean_optical_depth]
        )
        
        interpolator = interp1d(
            self.tau_integ[self.simulation_state.geometry.v_inner_boundary_index:],
            self.simulation_state.geometry.v_inner_active,  # Only use the active values as we only need a numerical estimate, not an index
            fill_value="extrapolate",
        )

        estimated_v_inner = interpolator(self.TAU_TARGET) * u.cm / u.s

        if estimated_v_inner < self.simulation_state.geometry.v_inner[0]:
            estimated_v_inner = self.simulation_state.geometry.v_inner[0]
            logger.warning(
                "WARNING: v_inner_boundary outside of simulation, setting to first shell"
            )
        elif estimated_v_inner > self.simulation_state.geometry.v_inner[-1]:
            estimated_v_inner = self.simulation_state.geometry.v_inner[-1]
            logger.warning(
                "WARNING: v_inner_boundary outside of simulation, setting to last shell"
            )

        return estimated_v_inner

    @property
    def property_mask(self):
        mask = np.zeros(
            (len(self.simulation_state.geometry.r_inner)), dtype=bool
        )
        mask[
            self.simulation_state.geometry.v_inner_boundary_index : self.simulation_state.geometry.v_outer_boundary_index
        ] = True
        return mask

    def get_convergence_estimates(self, transport_state):
        """Compute convergence estimates from the transport state

        Parameters
        ----------
        transport_state : MonteCarloTransportState
            Transport state object to compute estimates

        Returns
        -------
        dict
            Convergence estimates
        EstimatedRadiationFieldProperties
            Dilute radiation file and j_blues dataclass
        """
        estimates = super().get_convergence_estimates(transport_state)

        estimated_v_inner = self.estimate_v_inner()

        estimates[0].update(
            {"v_inner_boundary": estimated_v_inner, "mask": self.property_mask}
        )

        return estimates

    def reproject(self, array_one, mask_one, array_two, mask_two):
        """Reprojects two sub_arrays defined by a set of masks onto an array where the masks of both objects are true

        let A1, A2 be arrays of size gemetry.no_of_shells and
            a1 = A1[mask_one],
            a2 = A2[mask_two]
        find a1*, a2* s.t.
            a1* = A1[mask_one & mask_two],
            a2* = A2[mask_one & mask_two]
        this is equivalent to
            a1* = A1[mask_one][mask_two[mask_one]] = a1[mask_two[mask_one]],
            a2* = A2[mask_two][mask_one[mask_two]] = a2[mask_one[mask_two]]

        Parameters
        ----------
        array1 : np.ndarray
            A sub array of an array with the shape of a geometry property
        mask_one : np.ndarray(bool)
            Mask such that the parrent array accessed at this mask gives a1
        array_two : np.ndarray
            A sub array of an array with the shape of a geometry property
        mask_two : np.ndarray(bool)
            Mask such that the parrent array accessed at this mask gives a2

        Returns
        -------
        array_one*
            reprojection of array_one onto mask_one & mask_two
        array_two*
            reprojection of array_two onto mask_one & mask_two
        """
        return array_one[mask_two[mask_one]], array_two[mask_one[mask_two]]

    def print_mask(self, mask):
        return "".join([{True: "-", False: "X"}[m] for m in mask]).join("[]")

    def check_convergence(
        self,
        estimated_values,
    ):
        """Check convergence status for a dict of estimated values

        Parameters
        ----------
        estimated_values : dict
            Estimates to check convergence

        Returns
        -------
        bool
            If convergence has occurred
        """
        convergence_statuses = []

        mask = estimated_values["mask"]

        if not np.all(mask == self.property_mask):
            convergence_statuses.append(False)
            # Need to set status to False if change in mask size
            logger.info(
                f"Resized Geometry, Convergence Suppressed\n"
                f"\t  Old Geometry: {self.print_mask(mask)}\n"
                f"\t  New Geometry: {self.print_mask(self.property_mask)}"
            )

        for key, solver in self.convergence_solvers.items():
            current_value = getattr(self.simulation_state, key)
            estimated_value = estimated_values[key]

            if key in ["t_radiative", "dilution_factor"]:
                current_value, estimated_value = self.reproject(
                    current_value, self.property_mask, estimated_value, mask
                )
                no_of_shells = current_value.shape[0]
            else:
                no_of_shells = 1

            convergence_statuses.append(
                solver.get_convergence_status(
                    current_value, estimated_value, no_of_shells
                )
            )

        if np.all(convergence_statuses):
            hold_iterations = self.convergence_strategy.hold_iterations
            self.consecutive_converges_count += 1
            logger.info(
                f"Iteration converged {self.consecutive_converges_count:d}/{(hold_iterations + 1):d} consecutive "
                f"times."
            )

            return self.consecutive_converges_count == hold_iterations + 1

        self.consecutive_converges_count = 0
        return False

    def solve_simulation_state(self, estimated_values):
        """Update the simulation state with new inputs computed from previous
        iteration estimates.

        Parameters
        ----------
        estimated_values : dict
            Estimated from the previous iterations

        Returns
        -------
        next_values : dict
            The next values assigned to the simulation state
        """
        next_values = super().solve_simulation_state(estimated_values)
        self.simulation_state.geometry.v_inner_boundary = next_values[
            "v_inner_boundary"
        ]
        self.simulation_state.blackbody_packet_source.radius = (
            self.simulation_state.r_inner[0]
        )

        return next_values

    def solve_plasma(
        self,
        estimated_radfield_properties,
        mask,
    ):
        """Update the plasma solution with the new radiation field estimates

        Parameters
        ----------
        estimated_radfield_properties : EstimatedRadiationFieldProperties
            The radiation field properties to use for updating the plasma
        radiation_field: tardis.plasma.radiation_field.RadiationField
            Current radiation field object from the last iteration

        Raises
        ------
        ValueError
            If the plasma solver radiative rates type is unknown
        """
        radiation_field = DilutePlanckianRadiationField(
            temperature=self.simulation_state.radiation_field_state.temperature,
            dilution_factor=self.simulation_state.radiation_field_state.dilution_factor,
        )
        update_properties = dict(
            dilute_planckian_radiation_field=radiation_field
        )
        # A check to see if the plasma is set with JBluesDetailed, in which
        # case it needs some extra kwargs.
        if (
            self.plasma_solver.plasma_solver_settings.RADIATIVE_RATES_TYPE
            == "blackbody"
        ):
            planckian_radiation_field = (
                radiation_field.to_planckian_radiation_field()
            )
            j_blues = planckian_radiation_field.calculate_mean_intensity(
                self.plasma_solver.atomic_data.lines.nu.values
            )
            update_properties["j_blues"] = pd.DataFrame(
                j_blues, index=self.plasma_solver.atomic_data.lines.index
            )
        elif (
            self.plasma_solver.plasma_solver_settings.RADIATIVE_RATES_TYPE
            == "dilute-blackbody"
        ):
            j_blues = radiation_field.calculate_mean_intensity(
                self.plasma_solver.atomic_data.lines.nu.values
            )
            update_properties["j_blues"] = pd.DataFrame(
                j_blues, index=self.plasma_solver.atomic_data.lines.index
            )
        elif (
            self.plasma_solver.plasma_solver_settings.RADIATIVE_RATES_TYPE
            == "detailed"
        ):
            j_blues = radiation_field.calculate_mean_intensity(
                self.plasma_solver.atomic_data.lines.nu.values
            )
            j_blues[:, mask] = estimated_radfield_properties.j_blues
            update_properties["j_blues"] = pd.DataFrame(
                j_blues,
                index=self.plasma_solver.atomic_data.lines.index,
            )
        else:
            raise ValueError(
                f"radiative_rates_type type unknown - {self.plasma.plasma_solver_settings.RADIATIVE_RATES_TYPE}"
            )

        self.plasma_solver.update(**update_properties)

    def run(self):
        """Run the TARDIS simulation until convergence is reached"""
        self.converged = False
        while self.completed_iterations < self.total_iterations - 1:
            logger.info(
                f"\n\tStarting iteration {(self.completed_iterations + 1):d} of {self.total_iterations:d}"
            )
            self.store_plasma_state(
                self.completed_iterations,
                self.simulation_state.t_radiative,
                self.simulation_state.dilution_factor,
                self.plasma_solver.electron_densities,
                self.simulation_state.t_inner,
                self.simulation_state.geometry.v_inner_boundary,
                self.tau_integ
            )

            # Note that we are updating the class attribute here to ensure consistency
            self.opacity_states = self.solve_opacity()

            transport_state, virtual_packet_energies = self.solve_montecarlo(
                self.opacity_states, self.real_packet_count
            )

            (
                estimated_values,
                estimated_radfield_properties,
            ) = self.get_convergence_estimates(transport_state)

            self.solve_simulation_state(estimated_values)

            self.solve_plasma(
                estimated_radfield_properties,
                estimated_values["mask"],
            )

            self.converged = self.check_convergence(estimated_values)

            self.completed_iterations += 1

            if self.converged and self.convergence_strategy.stop_if_converged:
                break

        if self.converged:
            logger.info("\n\tStarting final iteration")
        else:
            logger.error(
                "\n\tITERATIONS HAVE NOT CONVERGED, starting final iteration"
            )
        self.opacity_states = self.solve_opacity()
        transport_state, virtual_packet_energies = self.solve_montecarlo(
            self.opacity_states,
            self.final_iteration_packet_count,
            self.virtual_packet_count,
        )

        self.store_plasma_state(
                self.completed_iterations,
                self.simulation_state.t_radiative,
                self.simulation_state.dilution_factor,
                self.plasma_solver.electron_densities,
                self.simulation_state.t_inner,
                self.simulation_state.geometry.v_inner_boundary,
                self.tau_integ
            )
        
        self.reshape_store_plasma_state(self.completed_iterations)
        
        self.initialize_spectrum_solver(
            transport_state,
            self.opacity_states,
            virtual_packet_energies,
        )
