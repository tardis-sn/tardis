import logging

import numpy as np
from astropy import units as u

from tardis.simulation.convergence import ConvergenceSolver
from tardis.workflows.simple_simulation import SimpleSimulation
from tardis.workflows.util import get_tau_integ
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from scipy.interpolate import interp1d
import pandas as pd

# logging support
logger = logging.getLogger(__name__)

# TODO:
# Option to estimate initial v_inner from electron opacity
# Add option for desired optical depth
# Think about csvy vs other formats for specifying grid
# Handle non-explicit formats when going out of the simulation


class InnerVelocitySimulationSolver(SimpleSimulation):

    TAU_TARGET = np.log(2.0 / 3)

    def __init__(self, configuration, mean_optical_depth="rossland", tau=None):

        super().__init__(configuration)
        self.mean_optical_depth = mean_optical_depth.lower()

        self.convergence_solvers["v_inner_boundary"] = ConvergenceSolver(
            self.convergence_strategy.v_inner_boundary
        )

        if tau is not None:
            self.TAU_TARGET = np.log(tau)

    def estimate_v_inner(self):
        """Compute the Rossland Mean Optical Depth,
        Estimate location where v_inner makes t=2/3 (or target)
        Extrapolate with exponential fits

        Need some way to return and inspect the optical depths for later logging"""
        pass

        tau_integ = np.log(
            get_tau_integ(
                self.plasma_solver,
                self.simulation_state,
            )[self.mean_optical_depth]
        )

        interpolator = interp1d(
            tau_integ,
            self.simulation_state.geometry.v_inner,  # Only use the active values as we only need a numerical estimate, not an index
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

    def reproject(self, a1, m1, a2, m2):
        """Reprojects two sub_arrays defined by a set of masks onto an array where the masks of both objects are true

        Parameters
        ----------
        a1 : np.ndarray
            A sub array of an array with the shape of a geometry property
        m1 : np.ndarray(bool)
            Mask such that the parrent array accessed at this mask gives a1
        a2 : np.ndarray
            A sub array of an array with the shape of a geometry property
        m2 : np.ndarray(bool)
            Mask such that the parrent array accessed at this mask gives a2

        Returns
        -------
        a1_joint
            reprojection of a1 onto m1 & m2
        a2_joint
        reprojection of a2 onto m1 & m2
        """

        a1_expanded = np.empty_like(
            a1,
            shape=(self.simulation_state.geometry.no_of_shells,),
        )
        a2_expanded = np.empty_like(a1_expanded)

        a1_expanded[m1] = a1
        a2_expanded[m2] = a2

        joint_mask = m1 & m2

        return a1_expanded[joint_mask], a2_expanded[joint_mask]

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
            logger.info(f"Resized Geometry, Convergence Suppressed")

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

    def solve_radiation_field(self):

        radiation_field = DilutePlanckianRadiationField(
            temperature=self.simulation_state.radiation_field_state.temperature,
            dilution_factor=self.simulation_state.radiation_field_state.dilution_factor,
        )

        return radiation_field

    def solve_plasma(
        self,
        estimated_radfield_properties,
        mask,
        radiation_field,
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
        converged = False
        while self.completed_iterations < self.total_iterations - 1:

            transport_state, virtual_packet_energies = self.solve_montecarlo(
                self.real_packet_count
            )

            (
                estimated_values,
                estimated_radfield_properties,
            ) = self.get_convergence_estimates(transport_state)

            self.solve_simulation_state(estimated_values)

            radiation_field = self.solve_radiation_field()

            self.solve_plasma(
                estimated_radfield_properties,
                estimated_values["mask"],
                radiation_field,
            )

            converged = self.check_convergence(estimated_values)

            self.completed_iterations += 1

            if converged and self.convergence_strategy.stop_if_converged:
                break

        if converged:
            logger.info("\n\tStarting final iteration")
        else:
            logger.error(
                "\n\tITERATIONS HAVE NOT CONVERGED, starting final iteration"
            )

        transport_state, virtual_packet_energies = self.solve_montecarlo(
            self.final_iteration_packet_count, self.virtual_packet_count
        )
        self.initialize_spectrum_solver(
            transport_state,
            virtual_packet_energies,
        )
