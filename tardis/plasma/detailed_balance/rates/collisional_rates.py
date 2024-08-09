import pandas as pd

from tardis.plasma.detailed_balance.rates.collision_strengths import (
    UpsilonCMFGENSolver,
    UpsilonRegemorterSolver,
)


class ThermalCollisionalRateSolver:
    def __init__(
        self,
        radiative_transitions,
        thermal_collisional_strength_temperatures,
        thermal_collisional_strengths,
        collision_strength_type,
        collisional_strength_approximation="regemorter",
    ):

        if collision_strength_type == "cmfgen":
            self.thermal_collision_strength_solver = UpsilonCMFGENSolver(
                thermal_collisional_strength_temperatures,
                thermal_collisional_strengths,
            )

        # find the transitions that have radiative rate data but no collisional data
        approximate_collisional_strength_index = (
            radiative_transitions.index.difference(
                thermal_collisional_strengths.index
            )
        )

        if collisional_strength_approximation == "regemorter":
            self.thermal_collision_strength_approximator = (
                UpsilonRegemorterSolver(
                    radiative_transitions.loc[
                        approximate_collisional_strength_index
                    ]
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

        return pd.concat(
            [
                thermal_collision_strengths,
                thermal_collision_strength_approximated,
            ]
        ).sort_index()
