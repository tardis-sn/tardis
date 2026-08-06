from collections.abc import Callable

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from tardis.plasma.exceptions import PlasmaIonizationError


class ChargeConservationSolver:
    """Coordinate independent elemental population solves by charge conservation.

    At a fixed trial electron density, the elemental population solves are
    independent. The coordinator combines their absolute ion populations and
    solves

    ``Q(n_e) = sum(Z * N_Z,j) - n_e = 0``

    independently for each shell.

    Parameters
    ----------
    elemental_number_density : pandas.DataFrame
        Elemental number densities indexed by atomic number and columned by
        shell.
    tolerance : float, optional
        Maximum allowed normalized charge residual, by default ``1e-10``.
    max_solver_iterations : int, optional
        Maximum number of iterations for each bounded scalar solve, by default
        ``100``.
    """

    def __init__(
        self,
        elemental_number_density: pd.DataFrame,
        tolerance: float = 1e-10,
        max_solver_iterations: int = 100,
    ) -> None:
        """Initialize the charge-conservation coordinator.

        Parameters
        ----------
        elemental_number_density : pandas.DataFrame
            Elemental number densities indexed by atomic number and columned
            by shell.
        tolerance : float, optional
            Maximum normalized charge residual accepted at convergence.
        max_solver_iterations : int, optional
            Maximum number of Brent iterations for one shell.
        """
        density = elemental_number_density.astype(float)
        self.elemental_number_density = density
        self.tolerance = tolerance
        self.max_solver_iterations = max_solver_iterations

    def _calculate_charge_residual(
        self,
        electron_number_density: pd.Series,
        solve_trial_populations: Callable[[pd.Series], pd.DataFrame],
    ) -> tuple[pd.DataFrame, pd.Series]:
        """Calculate the unnormalized charge residual from trial populations.

        Parameters
        ----------
        electron_number_density : pandas.Series
            Trial electron densities indexed by shell.
        solve_trial_populations : callable
            Fixed-electron-density population calculation used for this solve.

        Returns
        -------
        pandas.DataFrame
            Absolute ion populations returned by the fixed-density solve.
        pandas.Series
            Charge density minus electron density for each shell.

        Raises
        ------
        PlasmaIonizationError
            If the charge residual is nonfinite.
        """
        ion_populations = solve_trial_populations(
            electron_number_density.copy()
        )
        # Weight every retained ion stage by its charge, including the bare
        # nucleus, before subtracting the common electron density.
        ion_numbers = ion_populations.index.get_level_values("ion_number")
        charge_density = ion_populations.mul(ion_numbers, axis=0).sum(axis=0)
        charge_residual = charge_density - electron_number_density
        if not np.isfinite(charge_residual.to_numpy()).all():
            raise PlasmaIonizationError(
                "Charge-conservation residual is nonfinite"
            )
        return ion_populations, charge_residual

    def _calculate_shell_charge_residual(
        self,
        electron_number_density: float,
        shell: int,
        solve_trial_populations: Callable[[pd.Series], pd.DataFrame],
    ) -> float:
        """Calculate one shell's charge residual at a trial density.

        Parameters
        ----------
        electron_number_density : float
            Trial electron density for ``shell``.
        shell : int
            Shell whose electron density is being evaluated.
        solve_trial_populations : callable
            Fixed-electron-density population calculation used for this solve.

        Returns
        -------
        float
            Charge density minus electron density for ``shell``.
        """
        trial_density = pd.Series(
            0.0, index=self.elemental_number_density.columns
        )
        # Shell equations are independent, so other shells can remain at zero
        # while Brent's method evaluates this shell.
        trial_density[shell] = electron_number_density
        _, charge_residual = self._calculate_charge_residual(
            trial_density, solve_trial_populations
        )
        return charge_residual[shell]

    def _solve_shell(
        self,
        shell: int,
        maximum: float,
        solve_trial_populations: Callable[[pd.Series], pd.DataFrame],
    ) -> float:
        """Solve one shell on the physical electron-density interval.

        Parameters
        ----------
        shell : int
            Shell whose electron density is being solved.
        maximum : float
            Maximum physical electron density when every element is fully
            ionized in this shell.
        solve_trial_populations : callable
            Fixed-electron-density population calculation used for this solve.

        Returns
        -------
        float
            Charge-conserving electron density for ``shell``.

        Raises
        ------
        PlasmaIonizationError
            If no bracket exists on the physical interval.
        """
        if maximum == 0.0:
            return 0.0

        lower_residual = self._calculate_shell_charge_residual(
            0.0, shell, solve_trial_populations
        )
        upper_residual = self._calculate_shell_charge_residual(
            maximum, shell, solve_trial_populations
        )
        if lower_residual * upper_residual > 0.0:
            raise PlasmaIonizationError(
                f"Charge-conservation residual does not bracket shell {shell}: "
                f"Q(0)={lower_residual}, Q(ne_max)={upper_residual}"
            )

        if lower_residual == 0.0:
            return 0.0
        if upper_residual == 0.0:
            return maximum
        return brentq(
            self._calculate_shell_charge_residual,
            0.0,
            maximum,
            args=(shell, solve_trial_populations),
            xtol=1e-14,
            rtol=1e-14,
            maxiter=self.max_solver_iterations,
        )

    def solve(
        self,
        solve_trial_populations: Callable[[pd.Series], pd.DataFrame],
    ) -> tuple[pd.Series, pd.DataFrame]:
        """Solve charge conservation independently in every shell.

        Parameters
        ----------
        solve_trial_populations : callable
            Fixed-electron-density population calculation used for this solve.
            It receives a copy of the trial electron densities and returns
            absolute ion populations indexed by ``(atomic_number,
            ion_number)``.

        Returns
        -------
        tuple[pandas.Series, pandas.DataFrame]
            Final electron density and absolute ion populations, respectively.

        Raises
        ------
        PlasmaIonizationError
            If a shell cannot be bracketed or the final normalized residual
            exceeds ``tolerance``.
        """
        columns = self.elemental_number_density.columns
        # The fully ionized abundance sets the upper physical bound for n_e.
        atomic_numbers = self.elemental_number_density.index.to_numpy(
            dtype=float
        )
        maximum = self.elemental_number_density.mul(
            atomic_numbers, axis=0
        ).sum()
        electron_number_density = pd.Series(index=columns, dtype=float)
        for shell in columns:
            electron_number_density[shell] = self._solve_shell(
                shell, maximum[shell], solve_trial_populations
            )

        # Recompute populations once at the shared final electron densities so
        # the returned populations and residuals describe the same state.
        ion_populations, charge_residual = self._calculate_charge_residual(
            electron_number_density, solve_trial_populations
        )
        normalized_charge_residual = charge_residual / maximum.replace(0.0, 1.0)
        if (np.abs(normalized_charge_residual) > self.tolerance).any():
            raise PlasmaIonizationError(
                "Final charge-conservation residual exceeds tolerance: "
                f"{normalized_charge_residual.to_dict()}"
            )
        return electron_number_density, ion_populations
