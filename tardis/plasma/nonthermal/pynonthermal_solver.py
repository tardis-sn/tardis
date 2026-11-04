"""Translate shell plasma states to fixed rates from pynonthermal."""

import logging
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from time import perf_counter
from typing import Any

import numpy as np
import pandas as pd
from astropy import units as u

from tardis.plasma.equilibrium.inputs import FixedRateState

logger = logging.getLogger(__name__)

BOUND_BOUND_INDEX_NAMES = [
    "atomic_number",
    "ion_number",
    "ion_number_source",
    "ion_number_destination",
    "level_number_source",
    "level_number_destination",
]
IONIZATION_INDEX_NAMES = [
    "atomic_number",
    "ion_number",
    "ion_number_source",
    "ion_number_destination",
]

# Pynonthermal uses the ARTIS He I level ordering. These destinations are the
# corresponding 37-level TARDIS helium superlevels used by the Boyle model.
PYNONTHERMAL_TO_TARDIS_HE_I_DESTINATIONS = {
    **{level_number: ((level_number, 1.0),) for level_number in range(22)},
    22: ((23, 1.0),),
    23: ((22, 1.0),),
    24: ((24, 21 / 48), (26, 27 / 48)),
    25: ((25, 7 / 16), (27, 9 / 16)),
    26: ((28, 1.0),),
    27: ((29, 1.0),),
    29: ((29, 1.0),),
    30: ((29, 1.0),),
    32: ((29, 1.0),),
    28: ((30, 1.0),),
    31: ((30, 1.0),),
    33: ((30, 1.0),),
    34: ((30, 1.0),),
    35: ((31, 1.0),),
    37: ((31, 1.0),),
    38: ((31, 1.0),),
    40: ((31, 1.0),),
    36: ((32, 1.0),),
    39: ((32, 1.0),),
    41: ((32, 1.0),),
    42: ((32, 1.0),),
    44: ((33, 1.0),),
    43: ((34, 1.0),),
    46: ((35, 1.0),),
    45: ((36, 1.0),),
}


@dataclass(frozen=True)
class NonthermalPlasmaInput:
    """Shell state required by nonthermal degradation.

    Parameters
    ----------
    deposition_rate_density : astropy.units.Quantity
        Deposited power per unit volume in shell order.
    electron_number_density : astropy.units.Quantity
        Free-electron number density in shell order.
    electron_temperature : astropy.units.Quantity
        Electron temperature in shell order.
    ion_number_density : pandas.DataFrame
        Ion populations indexed by zero-based TARDIS ion number, with shells
        in columns.
    """

    deposition_rate_density: u.Quantity
    electron_number_density: u.Quantity
    electron_temperature: u.Quantity
    ion_number_density: pd.DataFrame


@dataclass(frozen=True)
class NonthermalPlasmaState:
    """Fixed TARDIS rates and energy fractions from nonthermal degradation."""

    fixed_rates: FixedRateState
    deposited_energy_fractions: pd.DataFrame


class PynonthermalHeliumSolver:
    """Calculate helium fixed rates from independent shell solves."""

    def __init__(
        self,
        destination_mapping: Mapping[
            int, tuple[tuple[int, float], ...]
        ] = PYNONTHERMAL_TO_TARDIS_HE_I_DESTINATIONS,
        solver_factory: Callable[..., Any] | None = None,
        minimum_energy: float = 1.0,
        maximum_energy: float = 3000.0,
        number_of_energy_points: int = 4000,
        minimum_ion_number_fraction: float = 0.0,
        energy_closure_tolerance: float = 5.0e-2,
    ) -> None:
        """Initialize the pynonthermal adapter.

        Parameters
        ----------
        destination_mapping : mapping
            Pynonthermal upper levels mapped to TARDIS destinations and
            branching fractions.
        solver_factory : callable, optional
            Spencer-Fano solver constructor. Primarily supports deterministic
            tests without importing the optional dependency.
        minimum_energy : float, optional
            Lower Spencer-Fano energy-grid boundary in eV.
        maximum_energy : float, optional
            Upper Spencer-Fano energy-grid boundary in eV.
        number_of_energy_points : int, optional
            Number of Spencer-Fano energy-grid points.
        minimum_ion_number_fraction : float, optional
            Minimum shell ion number fraction included as a competing
            degradation channel. Positive He I and He II populations are
            always retained because this adapter returns their rates. Zero
            preserves the complete composition.
        energy_closure_tolerance : float, optional
            Allowed absolute departure of reported energy fractions from one.
        """
        self.destination_mapping = dict(destination_mapping)
        self.solver_factory = solver_factory
        self.minimum_energy = minimum_energy
        self.maximum_energy = maximum_energy
        self.number_of_energy_points = number_of_energy_points
        if not 0.0 <= minimum_ion_number_fraction <= 1.0:
            raise ValueError(
                "minimum_ion_number_fraction must be between zero and one."
            )
        self.minimum_ion_number_fraction = minimum_ion_number_fraction
        self.energy_closure_tolerance = energy_closure_tolerance
        for destinations in self.destination_mapping.values():
            branching_fractions = np.asarray(
                [fraction for _, fraction in destinations], dtype=float
            )
            if np.any(branching_fractions < 0.0) or not np.isclose(
                branching_fractions.sum(), 1.0
            ):
                raise ValueError(
                    "Every pynonthermal destination mapping must have "
                    "nonnegative branching fractions that sum to one."
                )

    def _get_solver_factory(self) -> Callable[..., Any]:
        """Return the configured or lazily imported solver constructor."""
        if self.solver_factory is not None:
            return self.solver_factory
        try:
            from pynonthermal import SpencerFanoSolver  # noqa: PLC0415
        except ImportError as exc:
            raise ImportError(
                "PynonthermalHeliumSolver requires the optional "
                "'pynonthermal' package."
            ) from exc
        return SpencerFanoSolver

    @staticmethod
    def _validate_input(plasma_input: NonthermalPlasmaInput) -> None:
        """Validate shell alignment and physical input domains."""
        shell_count = len(plasma_input.ion_number_density.columns)
        arrays = {
            "deposition_rate_density": plasma_input.deposition_rate_density.to_value(
                u.erg / u.s / u.cm**3
            ),
            "electron_number_density": plasma_input.electron_number_density.to_value(
                u.cm**-3
            ),
            "electron_temperature": plasma_input.electron_temperature.to_value(
                u.K
            ),
        }
        for name, values in arrays.items():
            values = np.atleast_1d(values)
            if values.shape != (shell_count,):
                raise ValueError(f"{name} must contain one value per shell.")
            if not np.all(np.isfinite(values)) or np.any(values < 0.0):
                raise ValueError(f"{name} must be finite and nonnegative.")
        positive_deposition = arrays["deposition_rate_density"] > 0.0
        if np.any(
            positive_deposition & (arrays["electron_number_density"] <= 0.0)
        ):
            raise ValueError(
                "electron_number_density must be positive in depositing shells."
            )
        ion_values = plasma_input.ion_number_density.to_numpy(dtype=float)
        if not np.all(np.isfinite(ion_values)) or np.any(ion_values < 0.0):
            raise ValueError(
                "ion_number_density must be finite and nonnegative."
            )

    def solve(
        self, plasma_input: NonthermalPlasmaInput
    ) -> NonthermalPlasmaState:
        """Solve every positive-deposition shell and return fixed rates."""
        self._validate_input(plasma_input)
        shell_columns = plasma_input.ion_number_density.columns
        deposition_ev = np.atleast_1d(
            plasma_input.deposition_rate_density.to_value(u.eV / u.s / u.cm**3)
        )
        deposition_erg = np.atleast_1d(
            plasma_input.deposition_rate_density.to_value(u.erg / u.s / u.cm**3)
        )
        electron_density = np.atleast_1d(
            plasma_input.electron_number_density.to_value(u.cm**-3)
        )
        electron_temperature = np.atleast_1d(
            plasma_input.electron_temperature.to_value(u.K)
        )
        excitation_rates: dict[tuple[int, ...], np.ndarray] = {}
        ionization_index = pd.MultiIndex.from_tuples(
            [(2, 0, 0, 1), (2, 1, 1, 2)],
            names=IONIZATION_INDEX_NAMES,
        )
        ionization_rates = pd.DataFrame(
            0.0, index=ionization_index, columns=shell_columns
        )
        heating_rate_density = pd.Series(
            0.0, index=shell_columns, name="nonthermal_heating_rate_density"
        )
        fractions = pd.DataFrame(
            0.0,
            index=["heating", "excitation", "ionization"],
            columns=shell_columns,
        )

        for shell_idx, shell in enumerate(shell_columns):
            if deposition_ev[shell_idx] == 0.0:
                logger.info(
                    "Skipping pynonthermal shell %s (%d/%d): zero deposition",
                    shell,
                    shell_idx + 1,
                    len(shell_columns),
                )
                continue
            shell_ion_population = plasma_input.ion_number_density[shell]
            total_ion_number_density = float(shell_ion_population.sum())
            eligible_ions = []
            eligible_ion_number_density = 0.0
            retained_ion_number_density = 0.0
            candidate_ion_count = 0
            for (
                atomic_number,
                ion_number,
            ), number_density_value in shell_ion_population.items():
                number_density = float(number_density_value)
                if number_density == 0.0 or ion_number >= atomic_number:
                    continue
                candidate_ion_count += 1
                eligible_ion_number_density += number_density
                required_helium_ion = atomic_number == 2 and ion_number in (
                    0,
                    1,
                )
                ion_number_fraction = (
                    number_density / total_ion_number_density
                )
                if (
                    not required_helium_ion
                    and ion_number_fraction
                    < self.minimum_ion_number_fraction
                ):
                    continue
                eligible_ions.append(
                    (int(atomic_number), int(ion_number), number_density)
                )
                retained_ion_number_density += number_density
            shell_start_time = perf_counter()
            logger.info(
                "Starting pynonthermal shell %s (%d/%d) with %d/%d ion "
                "stages (retained eligible-ion density fraction %.8f)",
                shell,
                shell_idx + 1,
                len(shell_columns),
                len(eligible_ions),
                candidate_ion_count,
                (
                    retained_ion_number_density / eligible_ion_number_density
                    if eligible_ion_number_density > 0.0
                    else 0.0
                ),
            )
            solver = self._get_solver_factory()(
                emin_ev=self.minimum_energy,
                emax_ev=self.maximum_energy,
                npts=self.number_of_energy_points,
                verbose=False,
            )
            for atomic_number, ion_number, number_density in eligible_ions:
                ion_stage = int(ion_number) + 1
                solver.add_ionisation(
                    Z=int(atomic_number),
                    ion_stage=ion_stage,
                    n_ion=number_density,
                )
                try:
                    solver.add_ion_ltepopexcitation(
                        Z=int(atomic_number),
                        ion_stage=ion_stage,
                        n_ion=number_density,
                        temperature=float(electron_temperature[shell_idx]),
                    )
                except AssertionError as exc:
                    # Pynonthermal's bundled excitation database does not
                    # contain every ion supported by its ionization data.
                    if "No excitation data" not in str(exc):
                        raise

            solver.solve(
                depositionratedensity_ev=float(deposition_ev[shell_idx]),
                override_n_e=float(electron_density[shell_idx]),
            )
            solver.analyse_ntspectrum()
            heating_fraction = float(solver.get_frac_heating())
            if not np.isfinite(heating_fraction) or heating_fraction < 0.0:
                raise FloatingPointError(
                    f"Invalid nonthermal heating fraction in shell {shell}."
                )
            heating_rate_density.loc[shell] = (
                heating_fraction * deposition_erg[shell_idx]
            )
            fractions.loc["heating", shell] = heating_fraction

            for ion_number in (0, 1):
                ion_stage = ion_number + 1
                if (2, ion_stage) not in solver.ionpopdict:
                    continue
                ionization_rate = float(
                    solver.get_ionisation_ratecoeff(2, ion_stage)
                )
                if not np.isfinite(ionization_rate) or ionization_rate < 0.0:
                    raise FloatingPointError(
                        f"Invalid helium ionization rate in shell {shell}."
                    )
                ionization_rates.loc[
                    (2, ion_number, ion_number, ion_number + 1), shell
                ] = ionization_rate

            for transition_key in solver.excitationlists.get((2, 1), {}):
                lower_level, pynonthermal_upper_level = transition_key
                destinations = self.destination_mapping.get(
                    pynonthermal_upper_level
                )
                if destinations is None:
                    continue
                rate_per_deposited_energy = float(
                    solver.get_excitation_ratecoeff(
                        Z=2,
                        ion_stage=1,
                        transitionkey=transition_key,
                    )
                )
                excitation_rate = (
                    rate_per_deposited_energy * deposition_ev[shell_idx]
                )
                if not np.isfinite(excitation_rate) or excitation_rate < 0.0:
                    raise FloatingPointError(
                        f"Invalid helium excitation rate in shell {shell}."
                    )
                for destination_level, branching_fraction in destinations:
                    transition = (
                        2,
                        0,
                        0,
                        0,
                        int(lower_level),
                        destination_level,
                    )
                    excitation_rates.setdefault(
                        transition, np.zeros(len(shell_columns))
                    )[shell_idx] += excitation_rate * branching_fraction

            if hasattr(solver, "get_frac_excitation_tot"):
                fractions.loc["excitation", shell] = float(
                    solver.get_frac_excitation_tot()
                )
            if hasattr(solver, "get_frac_ionisation_tot"):
                fractions.loc["ionization", shell] = float(
                    solver.get_frac_ionisation_tot()
                )
            if hasattr(solver, "get_frac_sum"):
                fraction_sum = float(solver.get_frac_sum())
                if (
                    not np.isfinite(fraction_sum)
                    or abs(fraction_sum - 1.0) > self.energy_closure_tolerance
                ):
                    raise FloatingPointError(
                        "Pynonthermal energy fractions do not close in shell "
                        f"{shell}: {fraction_sum:.6f}."
                    )
            logger.info(
                "Finished pynonthermal shell %s (%d/%d) in %.3f s",
                shell,
                shell_idx + 1,
                len(shell_columns),
                perf_counter() - shell_start_time,
            )

        bound_bound_index = pd.MultiIndex.from_tuples(
            excitation_rates, names=BOUND_BOUND_INDEX_NAMES
        )
        bound_bound_rates = pd.DataFrame(
            list(excitation_rates.values()),
            index=bound_bound_index,
            columns=shell_columns,
            dtype=float,
        )
        return NonthermalPlasmaState(
            FixedRateState(
                bound_bound_rates,
                ionization_rates,
                heating_rate_density,
            ),
            fractions,
        )
