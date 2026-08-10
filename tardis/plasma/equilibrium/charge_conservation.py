import numpy as np
import pandas as pd
from scipy.optimize import brentq

from tardis.plasma.equilibrium.population_state import PopulationState
from tardis.plasma.exceptions import PlasmaIonizationError


class ChargeConservationSolver:
    """Coordinate independent elemental population solves by charge conservation.

    At a chosen trial electron density, the elemental population solves are
    independent for each element. The coordinator combines their absolute
    ion populations and solves independently for each shell:

    ``Q(n_e) = sum_Z sum_j j * N_(Z,j)(n_e) - n_e = 0``

    Here, ``Q(n_e)`` is the charge-number residual in number-density units and
    ``n_e`` is the trial electron number density. The outer sum runs over
    elements identified by atomic number ``Z``; the inner sum runs over ion
    charges ``j`` from zero through ``Z``. ``N_(Z,j)(n_e)`` is the absolute
    number density of ion stage ``j`` for element ``Z``. Equivalently,
    ``N_(Z,j)(n_e) = N_Z * Y_(Z,j)(n_e)``, where ``N_Z`` is the total elemental
    number density and ``Y_(Z,j)(n_e)`` is the element-normalized ion-stage
    fraction.

    For numerical conditioning, each scalar solve uses
    ``x = n_e / n_e,max``, where
    ``n_e,max = sum_Z Z * N_Z`` is the fully ionized electron density, and
    ``Q_hat(x) = Q(x * n_e,max) / n_e,max`` on ``[0, 1]``.

    Parameters
    ----------
    elemental_number_density : pandas.DataFrame
        Elemental number densities indexed by atomic number and columned by
        shell.
    population_solver : object
        Solver used for every trial electron-density population evaluation.
    tolerance : float, optional
        Maximum allowed normalized charge residual, by default ``1e-10``.
    max_solver_iterations : int, optional
        Maximum number of iterations for each bounded scalar solve, by default
        ``100``.
    """

    def __init__(
        self,
        elemental_number_density: pd.DataFrame,
        population_solver: object,
        tolerance: float = 1e-10,
        max_solver_iterations: int = 100,
    ) -> None:
        """Initialize the charge-conservation solver.

        Parameters
        ----------
        elemental_number_density : pandas.DataFrame
            Elemental number densities indexed by atomic number and columned
            by shell.
        population_solver : object
            Solver exposing ``solve_charge_state_at_electron_density``.
        tolerance : float, optional
            Maximum normalized charge residual accepted at convergence.
        max_solver_iterations : int, optional
            Maximum number of Brent iterations for one shell.
        """
        density = elemental_number_density.astype(float)
        self.elemental_number_density = density
        self.population_solver = population_solver
        self.tolerance = tolerance
        self.max_solver_iterations = max_solver_iterations

    def _calculate_charge_residual(
        self,
        electron_number_density: pd.Series,
    ) -> tuple[PopulationState, pd.Series]:
        """Calculate the unnormalized charge residual from trial populations.

        Parameters
        ----------
        electron_number_density : pandas.Series
            Trial electron densities indexed by shell.

        Returns
        -------
        PopulationState
            Complete state returned for the trial electron density.
        pandas.Series
            Charge density minus electron density for each shell.

        Raises
        ------
        PlasmaIonizationError
            If the charge residual is nonfinite.
        """
        population_state = (
            self.population_solver.solve_charge_state_at_electron_density(
                electron_number_density.copy()
            )
        )
        ion_populations = population_state.ion_number_density
        # Weight every retained ion stage by its charge, including the bare
        # nucleus, before subtracting the common electron density.
        ion_numbers = ion_populations.index.get_level_values("ion_number")
        charge_density = ion_populations.mul(ion_numbers, axis=0).sum(axis=0)
        charge_residual = charge_density - electron_number_density
        if not np.isfinite(charge_residual.to_numpy()).all():
            raise PlasmaIonizationError(
                "Charge-conservation residual is nonfinite"
            )
        return population_state, charge_residual

    def _calculate_shell_charge_residual(
        self,
        electron_density_fraction: float,
        shell: int,
        maximum_electron_number_density: float,
    ) -> float:
        """Calculate one shell's normalized charge residual.

        Parameters
        ----------
        electron_density_fraction : float
            Trial electron density divided by its fully ionized upper bound.
        shell : int
            Shell whose electron density is being evaluated.
        maximum_electron_number_density : float
            Fully ionized electron density for ``shell``.

        Returns
        -------
        float
            Charge density minus electron density, divided by the fully
            ionized electron density for ``shell``.
        """
        trial_density = pd.Series(
            [electron_density_fraction * maximum_electron_number_density],
            index=pd.Index(
                [shell], name=self.elemental_number_density.columns.name
            ),
        )
        _, charge_residual = self._calculate_charge_residual(
            trial_density,
        )
        return charge_residual[shell] / maximum_electron_number_density

    def _solve_shell(
        self,
        shell: int,
        maximum_electron_number_density: float,
        electron_number_density_seed: float | None = None,
    ) -> float:
        """Solve one shell on the normalized electron-density interval.

        Parameters
        ----------
        shell : int
            Shell whose electron density is being solved.
        maximum_electron_number_density : float
            Maximum physical electron density when every element is fully
            ionized in this shell.
        electron_number_density_seed : float, optional
            Previous electron density used to choose the first sub-bracket.

        Returns
        -------
        float
            Charge-conserving electron density for ``shell``.

        Raises
        ------
        PlasmaIonizationError
            If no bracket exists on the physical interval.
        """
        if maximum_electron_number_density == 0.0:
            return 0.0

        lower_fraction = 0.0
        upper_fraction = 1.0
        lower_residual = None
        upper_residual = None
        if electron_number_density_seed is not None:
            seed_fraction = float(
                np.clip(
                    electron_number_density_seed
                    / maximum_electron_number_density,
                    0.0,
                    1.0,
                )
            )
            seed_residual = self._calculate_shell_charge_residual(
                seed_fraction,
                shell,
                maximum_electron_number_density,
            )
            if seed_residual == 0.0:
                return seed_fraction * maximum_electron_number_density
            if seed_residual > 0.0:
                lower_fraction = seed_fraction
                lower_residual = seed_residual
            else:
                upper_fraction = seed_fraction
                upper_residual = seed_residual

        if lower_residual is None:
            lower_residual = self._calculate_shell_charge_residual(
                lower_fraction,
                shell,
                maximum_electron_number_density,
            )
        if upper_residual is None:
            upper_residual = self._calculate_shell_charge_residual(
                upper_fraction,
                shell,
                maximum_electron_number_density,
            )
        if lower_residual * upper_residual > 0.0 and (
            lower_fraction != 0.0 or upper_fraction != 1.0
        ):
            lower_fraction = 0.0
            upper_fraction = 1.0
            lower_residual = self._calculate_shell_charge_residual(
                lower_fraction,
                shell,
                maximum_electron_number_density,
            )
            upper_residual = self._calculate_shell_charge_residual(
                upper_fraction,
                shell,
                maximum_electron_number_density,
            )
        if lower_residual * upper_residual > 0.0:
            raise PlasmaIonizationError(
                f"Charge-conservation residual does not bracket shell {shell}: "
                f"Q_hat(0)={lower_residual}, Q_hat(1)={upper_residual}"
            )

        if lower_residual == 0.0:
            return float(lower_fraction) * maximum_electron_number_density
        if upper_residual == 0.0:
            return upper_fraction * maximum_electron_number_density
        electron_density_fraction = brentq(
            self._calculate_shell_charge_residual,
            float(lower_fraction),
            float(upper_fraction),
            args=(
                shell,
                maximum_electron_number_density,
            ),
            xtol=1e-14,
            rtol=1e-14,
            maxiter=self.max_solver_iterations,
        )
        return electron_density_fraction * maximum_electron_number_density

    def solve(
        self,
        beta_sobolev: pd.DataFrame | None = None,
        electron_number_density_seed: pd.Series | None = None,
    ) -> PopulationState:
        """Solve charge conservation independently in every shell.

        Parameters
        ----------
        beta_sobolev : pandas.DataFrame, optional
            Sobolev escape probabilities indexed by line and columned by
            shell. Required only when the population solver prepares
            beta-dependent rates.
        electron_number_density_seed : pandas.Series, optional
            Previous electron densities used only to choose nearby Brent
            brackets. The full physical interval remains the fallback.

        Returns
        -------
        PopulationState
            Complete state evaluated once at the final charge-conserving
            electron densities.

        Raises
        ------
        PlasmaIonizationError
            If a shell cannot be bracketed or the final normalized residual
            exceeds ``tolerance``.
        """
        columns = self.elemental_number_density.columns
        prepare_solve = getattr(
            self.population_solver,
            "prepare_charge_conservation_solve",
            None,
        )
        if prepare_solve is not None:
            prepared_beta = (
                None if beta_sobolev is None else beta_sobolev.copy(deep=True)
            )
            prepare_solve(prepared_beta)
        # The fully ionized abundance sets the upper physical bound for n_e.
        atomic_numbers = self.elemental_number_density.index.to_numpy(
            dtype=float
        )
        maximum_electron_number_density = self.elemental_number_density.mul(
            atomic_numbers, axis=0
        ).sum()
        electron_number_density = pd.Series(index=columns, dtype=float)
        if electron_number_density_seed is not None:
            electron_number_density_seed = electron_number_density_seed.reindex(
                columns
            )
        for shell in columns:
            electron_number_density[shell] = self._solve_shell(
                shell,
                maximum_electron_number_density[shell],
                None
                if electron_number_density_seed is None
                else electron_number_density_seed[shell],
            )

        # Back-substitute once at the final roots so every returned population
        # and the validated residual describe the same state.
        population_state, charge_residual = self._calculate_charge_residual(
            electron_number_density
        )
        normalized_charge_residual = charge_residual / (
            maximum_electron_number_density.replace(0.0, 1.0)
        )
        if (np.abs(normalized_charge_residual) > self.tolerance).any():
            raise PlasmaIonizationError(
                "Final charge-conservation residual exceeds tolerance: "
                f"{normalized_charge_residual.to_dict()}"
            )
        return population_state
