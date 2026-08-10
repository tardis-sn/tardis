import numpy as np
import numpy.typing as npt
import pandas as pd
from astropy import units as u

from tardis.opacities.sobolev import (
    calculate_beta_sobolev,
    calculate_sobolev_line_opacity,
)
from tardis.plasma.equilibrium.charge_conservation import (
    ChargeConservationSolver,
)
from tardis.plasma.equilibrium.ion_populations import PopulationState
from tardis.plasma.exceptions import PlasmaIonizationError
from tardis.plasma.properties.radiative_properties import (
    calculate_stimulated_emission_factor,
)


class SobolevPopulationSolver:
    """Solve populations while refreshing population-dependent Sobolev rates.

    Parameters
    ----------
    charge_conservation_solver : ChargeConservationSolver
        Charge solver returning a complete coupled population state for a
        DataFrame of escape probabilities.
    lines : pandas.DataFrame
        Atomic line data indexed like the Sobolev property.
    time_explosion : astropy.units.Quantity
        Time since explosion used in the Sobolev optical depth.
    g : pandas.Series
        Statistical weights for all levels.
    lines_lower_level_index, lines_upper_level_index : array-like
        Positional indexes of the lower and upper levels for each line.
    metastability : pandas.Series
        Metastability flags for all levels.
    nlte_species : set[tuple[int, int]], optional
        Species for which population inversions are clipped to zero.
    tolerance : float, optional
        Relative beta and population update tolerance.
    max_iterations : int, optional
        Maximum number of population/beta updates.
    relaxation_floor : float, optional
        Smallest relaxation factor used after a growing update.
    """

    def __init__(
        self,
        charge_conservation_solver: ChargeConservationSolver,
        lines: pd.DataFrame,
        time_explosion: u.Quantity,
        g: pd.Series,
        lines_lower_level_index: npt.ArrayLike,
        lines_upper_level_index: npt.ArrayLike,
        metastability: pd.Series,
        nlte_species: set[tuple[int, int]] | None = None,
        tolerance: float = 1e-8,
        max_iterations: int = 100,
        relaxation_floor: float = 1e-3,
    ) -> None:
        self.charge_conservation_solver = charge_conservation_solver
        self.lines = lines
        self.time_explosion = time_explosion
        self.g = g
        self.lines_lower_level_index = lines_lower_level_index
        self.lines_upper_level_index = lines_upper_level_index
        self.metastability = metastability
        self.nlte_species = nlte_species
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.relaxation_floor = relaxation_floor

    @staticmethod
    def _relative_norm(previous: pd.DataFrame, current: pd.DataFrame) -> float:
        """Return a finite relative maximum norm for two aligned frames."""
        difference = np.abs(current.to_numpy() - previous.to_numpy())
        scale = np.maximum(np.abs(previous.to_numpy()), 1e-300)
        return float(np.max(difference / scale))

    def _calculate_sobolev_state(
        self, level_number_density: pd.DataFrame
    ) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
        """Calculate stimulated emission, optical depth, and beta."""
        stimulated_emission_factor = calculate_stimulated_emission_factor(
            self.g,
            level_number_density,
            self.lines_lower_level_index,
            self.lines_upper_level_index,
            self.metastability,
            self.lines,
            self.nlte_species,
        )
        tau_sobolevs = calculate_sobolev_line_opacity(
            self.lines,
            level_number_density,
            self.time_explosion,
            stimulated_emission_factor,
        )
        return (
            stimulated_emission_factor,
            tau_sobolevs,
            calculate_beta_sobolev(tau_sobolevs),
        )

    def solve(
        self,
        initial_level_number_density: pd.DataFrame,
        initial_beta_sobolev: pd.DataFrame | None = None,
    ) -> tuple[
        PopulationState,
        npt.NDArray[np.float64],
        pd.DataFrame,
        pd.DataFrame,
    ]:
        """Solve the coupled population and Sobolev escape-probability state.

        The supplied level populations are copied and used only as the initial
        state. For each iteration, the charge coordinator solves at the
        current escape probabilities. The resulting populations then
        determine the stimulated-emission factors, Sobolev optical depths, and
        proposed escape probabilities through plasma property calculations.
        The beta update is initially undamped. If its relative
        update norm grows, the relaxation factor is halved, down to
        ``relaxation_floor``, before the next population solve.

        Convergence requires both the relative beta update and the relative
        level-population update to be no larger than ``tolerance``. After
        convergence, the population solver is evaluated once more at the
        final escape probabilities. This back-substitution ensures that all
        returned values describe the same state and that the returned beta is
        reproducible from the returned optical depth.

        Parameters
        ----------
        initial_level_number_density : pandas.DataFrame
            Complete positive level state used to seed the iteration. Its
            index and columns define the level and shell layout returned by
            the population solver.
        initial_beta_sobolev : pandas.DataFrame, optional
            Escape probabilities from the previous iteration. If omitted,
            they are calculated from ``initial_level_number_density``. When
            provided, it must cover every line in ``lines`` and every shell in
            the initial level state.

        Returns
        -------
        tuple[PopulationState, numpy.ndarray, pandas.DataFrame, pandas.DataFrame]
            Values are returned in this order:

            * final coupled population state;
            * stimulated-emission factor for each line and shell;
            * Sobolev optical depth for each line and shell;
            * Sobolev escape probability for each line and shell.

        Raises
        ------
        PlasmaIonizationError
            If an update is nonfinite, the fixed point does not converge
            within ``max_iterations``, or final back-substitution is not
            self-consistent.
        """
        level_number_density = initial_level_number_density.copy(deep=True)
        if initial_beta_sobolev is None:
            _, _, beta_sobolev = self._calculate_sobolev_state(
                level_number_density
            )
        else:
            beta_sobolev = initial_beta_sobolev.reindex(
                self.lines.index, columns=level_number_density.columns
            ).copy()
        if not np.isfinite(beta_sobolev.to_numpy()).all():
            raise PlasmaIonizationError("Initial Sobolev beta is nonfinite")

        previous_update_norm = np.inf
        relaxation = 1.0
        converged = False
        for _iteration in range(1, self.max_iterations + 1):
            population_state = self.charge_conservation_solver.solve(
                beta_sobolev.copy(deep=True)
            )
            new_level_number_density = population_state.level_number_density
            (
                stimulated_emission_factor,
                tau_sobolevs,
                proposed_beta,
            ) = self._calculate_sobolev_state(new_level_number_density)
            update_norm = self._relative_norm(beta_sobolev, proposed_beta)
            population_norm = self._relative_norm(
                level_number_density, new_level_number_density
            )
            if not np.isfinite(update_norm + population_norm):
                raise PlasmaIonizationError(
                    "Sobolev fixed-point update is nonfinite"
                )
            if (
                update_norm <= self.tolerance
                and population_norm <= self.tolerance
            ):
                level_number_density = new_level_number_density
                beta_sobolev = proposed_beta
                converged = True
                break

            if update_norm > previous_update_norm:
                relaxation = max(self.relaxation_floor, relaxation / 2.0)
            next_beta = beta_sobolev + relaxation * (
                proposed_beta - beta_sobolev
            )
            level_number_density = new_level_number_density
            beta_sobolev = next_beta
            previous_update_norm = update_norm

        if not converged:
            raise PlasmaIonizationError(
                "Sobolev fixed-point iteration did not converge after "
                f"{self.max_iterations} iterations"
            )

        # Back-substitute once at the final beta so the returned population and
        # public Sobolev properties describe one identical state.
        population_state = self.charge_conservation_solver.solve(
            beta_sobolev.copy(deep=True)
        )
        level_number_density = population_state.level_number_density
        (
            stimulated_emission_factor,
            tau_sobolevs,
            final_beta,
        ) = self._calculate_sobolev_state(level_number_density)
        if self._relative_norm(beta_sobolev, final_beta) > self.tolerance:
            raise PlasmaIonizationError(
                "Final Sobolev back-substitution is not self-consistent"
            )
        return (
            population_state,
            stimulated_emission_factor,
            tau_sobolevs,
            final_beta,
        )
