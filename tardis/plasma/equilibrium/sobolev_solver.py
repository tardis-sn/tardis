import numpy as np
import numpy.typing as npt
import pandas as pd
from astropy import units as u

from tardis.opacities.sobolev import (
    calculate_beta_sobolev,
    calculate_sobolev_line_opacity,
    calculate_stimulated_emission_factor,
)
from tardis.plasma.equilibrium.charge_conservation import (
    ChargeConservationSolver,
)
from tardis.plasma.equilibrium.population_state import PopulationState
from tardis.plasma.exceptions import PlasmaIonizationError


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
        Relative beta and population update tolerance, by default
        ``1.49012e-8``, matching the legacy MINPACK level solve.
    max_iterations : int, optional
        Maximum number of population/beta updates, by default ``200``.
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
        tolerance: float = 1.49012e-8,
        max_iterations: int = 200,
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
        scale = np.maximum.reduce(
            (
                np.abs(previous.to_numpy()),
                np.abs(current.to_numpy()),
                np.full(previous.shape, 1e-300),
            )
        )
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
        initial_electron_densities: pd.Series | None = None,
    ) -> tuple[
        PopulationState,
        npt.NDArray[np.float64],
        pd.DataFrame,
        pd.DataFrame,
    ]:
        """Converge charge, populations, and Sobolev escape probabilities.

        Parameters
        ----------
        initial_level_number_density : pandas.DataFrame
            Complete previous level state used to seed the iteration.
        initial_beta_sobolev : pandas.DataFrame, optional
            Previous escape probabilities. When omitted, they are calculated
            from ``initial_level_number_density``.
        initial_electron_densities : pandas.Series, optional
            Previous electron densities used to seed nearby charge brackets.

        Returns
        -------
        tuple[PopulationState, numpy.ndarray, pandas.DataFrame, pandas.DataFrame]
            Converged populations, stimulated emission factors, Sobolev
            optical depths, and escape probabilities.

        Raises
        ------
        PlasmaIonizationError
            If the fixed point is nonfinite, fails to converge, or fails final
            back-substitution.
        """
        level_number_density = initial_level_number_density.copy(deep=True)
        if initial_beta_sobolev is None:
            _, _, beta_sobolev = self._calculate_sobolev_state(
                level_number_density
            )
        else:
            beta_sobolev = initial_beta_sobolev.reindex(
                self.lines.index, columns=level_number_density.columns
            ).copy(deep=True)
        if not np.isfinite(beta_sobolev.to_numpy()).all():
            raise PlasmaIonizationError("Initial Sobolev beta is nonfinite")

        accepted_merit = np.inf
        update_norm = np.inf
        population_norm = np.inf
        relaxation = 1.0
        retry_origin_beta = beta_sobolev
        retry_target_beta = beta_sobolev
        converged = False
        for _iteration in range(1, self.max_iterations + 1):
            population_state = self.charge_conservation_solver.solve(
                beta_sobolev.copy(deep=True),
                initial_electron_densities,
            )
            new_level_number_density = population_state.level_number_density
            _, _, proposed_beta = self._calculate_sobolev_state(
                new_level_number_density
            )
            update_norm = self._relative_norm(beta_sobolev, proposed_beta)
            population_norm = self._relative_norm(
                level_number_density, new_level_number_density
            )
            if not np.isfinite(update_norm + population_norm):
                raise PlasmaIonizationError(
                    "Sobolev fixed-point update is nonfinite"
                )
            update_merit = max(update_norm, population_norm)
            if update_merit > accepted_merit:
                if relaxation <= self.relaxation_floor:
                    raise PlasmaIonizationError(
                        "Sobolev fixed-point damping reached its floor: "
                        f"update_norm={update_norm}, "
                        f"population_norm={population_norm}, "
                        f"accepted_merit={accepted_merit}, "
                        f"relaxation={relaxation}"
                    )
                relaxation = max(self.relaxation_floor, relaxation / 2.0)
                beta_sobolev = retry_origin_beta + relaxation * (
                    retry_target_beta - retry_origin_beta
                )
                continue
            initial_electron_densities = (
                population_state.electron_densities.copy(deep=True)
            )
            if (
                update_norm <= self.tolerance
                and population_norm <= self.tolerance
            ):
                level_number_density = new_level_number_density
                beta_sobolev = proposed_beta
                converged = True
                break
            level_number_density = new_level_number_density
            retry_origin_beta = beta_sobolev
            retry_target_beta = proposed_beta
            beta_sobolev = proposed_beta
            accepted_merit = update_merit
            relaxation = 1.0

        if not converged:
            raise PlasmaIonizationError(
                "Sobolev fixed-point iteration did not converge after "
                f"{self.max_iterations} iterations: "
                f"update_norm={update_norm}, "
                f"population_norm={population_norm}, "
                f"relaxation={relaxation}"
            )

        population_state = self.charge_conservation_solver.solve(
            beta_sobolev.copy(deep=True),
            initial_electron_densities,
        )
        (
            stimulated_emission_factor,
            tau_sobolevs,
            final_beta,
        ) = self._calculate_sobolev_state(population_state.level_number_density)
        final_update_norm = self._relative_norm(beta_sobolev, final_beta)
        final_population_norm = self._relative_norm(
            level_number_density,
            population_state.level_number_density,
        )
        if (
            final_update_norm > self.tolerance
            or final_population_norm > self.tolerance
        ):
            raise PlasmaIonizationError(
                "Final Sobolev back-substitution is not self-consistent: "
                f"update_norm={final_update_norm}, "
                f"population_norm={final_population_norm}"
            )
        return (
            population_state,
            stimulated_emission_factor,
            tau_sobolevs,
            final_beta,
        )
