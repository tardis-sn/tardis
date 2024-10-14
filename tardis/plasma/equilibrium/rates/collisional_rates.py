import numpy as np
import pandas as pd
from astropy import units as u

from tardis import constants as const
from tardis.plasma.equilibrium.rates.collision_strengths import (
    UpsilonChiantiSolver,
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
        self.collision_strengths_type = collision_strengths_type
        if self.collision_strengths_type == "cmfgen":
            self.thermal_collision_strength_solver = UpsilonCMFGENSolver(
                thermal_collisional_strengths_temperatures,
                thermal_collisional_strengths,
            )
        elif self.collision_strengths_type == "chianti":
            self.thermal_collision_strength_solver = UpsilonChiantiSolver(
                thermal_collisional_strengths,
            )
        else:
            raise ValueError(
                f"collision_strengths_type {collision_strengths_type} not supported"
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

        self.g_l = self.levels.loc[
            self.all_collisional_strengths_index.droplevel("level_number_lower")
        ].g.values

        self.g_u = self.levels.loc[
            self.all_collisional_strengths_index.droplevel("level_number_upper")
        ].g.values

        if collisional_strength_approximation == "regemorter":
            self.thermal_collision_strength_approximator = (
                UpsilonRegemorterSolver(
                    radiative_transitions.loc[missing_collision_strengths_index]
                )
            )

    def solve(self, temperatures_electron):
        thermal_all_collision_strengths = self.calculate_collision_strengths(
            temperatures_electron
        )

        boltzmann_factor = np.exp(
            self.delta_energies[np.newaxis].T
            / (temperatures_electron * const.k_B),
        ).value
        collision_rates_coeff_lu = (
            (BETA_COLL / np.sqrt(temperatures_electron) * boltzmann_factor)
            .to("cm3 / s")
            .value
            * thermal_all_collision_strengths
        )  # see formula A2 in Przybilla, Butler 2004 - Apj 609, 1181

        collision_rates_coeff_ul = (
            (self.g_u / self.g_l)[np.newaxis].T
            / boltzmann_factor
            * collision_rates_coeff_lu
        )

        collision_rates_coeff_ul.index = (
            collision_rates_coeff_lu.index.swaplevel(
                "level_number_lower", "level_number_upper"
            )
        )

        collision_rates_coeff_df = pd.concat(
            [collision_rates_coeff_lu, collision_rates_coeff_ul]
        )
        collision_rates_coeff_df.ion_number_source = (
            collision_rates_coeff_df.ion_number
        )
        collision_rates_coeff_df.ion_number_destination = (
            collision_rates_coeff_df.ion_number
        )
        collision_rates_coeff_df.index.names = [
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
            "level_number_source",
            "level_number_destination",
        ]
        return collision_rates_coeff_df

    def calculate_collision_strengths(self, temperatures_electron):
        """
        Calculate collision strengths based on the provided electron temperatures.

        Parameters
        ----------
        temperatures_electron : array-like
            Array-like of electron temperatures.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the calculated collision strengths.
        """
        thermal_collision_strengths = (
            self.thermal_collision_strength_solver.solve(temperatures_electron)
        )
        thermal_collision_strength_approximated = (
            self.thermal_collision_strength_approximator.solve(
                temperatures_electron
            )
        )

        return pd.concat(
            [
                thermal_collision_strengths,
                thermal_collision_strength_approximated,
            ]
        ).sort_index()
