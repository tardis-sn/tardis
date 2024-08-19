from tardis.plasma.radiation_field.planck_rad_field import (
    DilutePlanckianRadiationField,
)
from tardis.simulation.convergence import ConvergenceSolver


class RadiationFieldConvergenceSolver:
    converge_separately = True
    temperature_strategy = {"damping_constant": 0, "threshold": 0, "type": "damped"}
    dilution_factor_strategy = {"damping_constant": 0, "threshold": 0, "type": "damped"}

    def __init__(self, strategy):
        """Create a RadiationFieldConvergenceSolver instance

        Parameters
        ----------
        strategy : dict
            Convergence strategy configuration including both t_rad and w for
            radiative temperature and dilution factor respectively.
        """
        self.temperature_strategy = strategy.t_rad
        self.dilution_factor_strategy = strategy.w

        self.temperature_solver = ConvergenceSolver(self.temperature_strategy)
        self.dilution_factor_solver = ConvergenceSolver(
            self.dilution_factor_strategy
        )

    def get_convergence_status(
        self, radiation_field, estimated_radiation_field, no_of_shells
    ):
        """Get the convergence status for the radiation field.

        Parameters
        ----------
        radiation_field : DilutePlanckianRadiationField
            The initial radiation field
        estimated_radiation_field : DilutePlanckianRadiationField
            The estimated radiation field
        no_of_shells : int
            Number of shells to average over

        Returns
        -------
        bool or Tuple of bool
            True if converged. Tuple is in order (temperature, dilution_factor)
        """
        temperature_converged = self.temperature_solver.get_convergence_status(
            radiation_field.temperature,
            estimated_radiation_field.temperature,
            no_of_shells,
        )
        dilution_factor_converged = (
            self.dilution_factor_solver.get_convergence_status(
                radiation_field.dilution_factor,
                estimated_radiation_field.dilution_factor,
                no_of_shells,
            )
        )

        if self.converge_separately:
            return temperature_converged, dilution_factor_converged

        return temperature_converged and dilution_factor_converged

    def converge(self, radiation_field, estimated_radiation_field):
        """Produce a new converged estimate for the radiation field.

        Parameters
        ----------
        radiation_field : DilutePlanckianRadiationField
            Initial radiation field
        estimated_radiation_field : DilutePlanckianRadiationField
            Estimated radiation field

        Returns
        -------
        DilutePlanckianRadiationField
            Converged radiation field
        """
        temperature_estimate = self.temperature_solver.converge(
            radiation_field.temperature, estimated_radiation_field.temperature
        )
        dilution_factor_estimate = self.dilution_factor_solver.converge(
            radiation_field.dilution_factor,
            estimated_radiation_field.dilution_factor,
        )

        return DilutePlanckianRadiationField(
            temperature_estimate, dilution_factor_estimate
        )
