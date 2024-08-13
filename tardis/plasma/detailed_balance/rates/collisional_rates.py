import numpy as np
import pandas as pd
from astropy import units as u

from tardis import constants as const
from tardis.plasma.detailed_balance.rates.collision_strengths import (
    UpsilonCMFGENSolver,
    UpsilonRegemorterSolver,
)

BETA_COLL = (
    (const.h**4 / (8 * const.k_B * const.m_e**3 * np.pi**3)) ** 0.5
).cgs


class ThermalCollisionalRateSolver:

    def __init__(
        self,
        levels,
        radiative_transitions,
        thermal_collisional_strengths_temperatures,
        thermal_collisional_strengths,
        collision_strengths_type,
        collisional_strength_approximation="regemorter",
    ):
        self.levels = levels
        if collision_strengths_type == "cmfgen":
            self.thermal_collision_strength_solver = UpsilonCMFGENSolver(
                thermal_collisional_strengths_temperatures,
                thermal_collisional_strengths,
            )
        self.radiative_transitions = radiative_transitions
        # find the transitions that have radiative rate data but no collisional data
        missing_collision_strengths_index = (
            radiative_transitions.index.difference(
                thermal_collisional_strengths.index
            )
        )
        self.all_collisional_strengths_index = (
            missing_collision_strengths_index.append(
                thermal_collisional_strengths.index
            ).sort_values()
        )
        self.delta_energies = (
            self.levels.loc[
                self.all_collisional_strengths_index.droplevel(
                    "level_number_upper"
                )
            ].energy.values
            - self.levels.loc[
                self.all_collisional_strengths_index.droplevel(
                    "level_number_lower"
                )
            ].energy.values
        ) * u.erg
        if collisional_strength_approximation == "regemorter":
            self.thermal_collision_strength_approximator = (
                UpsilonRegemorterSolver(
                    radiative_transitions.loc[missing_collision_strengths_index]
                )
            )

    def solve(self, temperatures_electron):
        thermal_collision_strengths = (
            self.thermal_collision_strength_solver.solve(temperatures_electron)
        )
        thermal_collision_strength_approximated = (
            self.thermal_collision_strength_approximator.solve(
                temperatures_electron
            )
        )

        thermal_all_collision_strengths = pd.concat(
            [
                thermal_collision_strengths,
                thermal_collision_strength_approximated,
            ]
        ).sort_index()

        boltzmann_factor = np.exp(
            self.delta_energies[np.newaxis].T
            / (temperatures_electron * const.k_B),
        ).value
        q_ij = (
            BETA_COLL / np.sqrt(temperatures_electron) * boltzmann_factor
        ).to(
            "cm3 / s"
        ).value * thermal_all_collision_strengths  # see formula A2 in Przybilla, Butler 2004 - Apj 609, 1181

        return q_ij
