from tardis.plasma.equilibrium.rates.collisional_ionization_strengths import (
    CollisionalIonizationSeaton,
)


class CollisionalIonizationRateSolver:
    """Solver for collisional ionization and recombination rates."""

    def __init__(self, photoionization_cross_sections):
        self.photoionization_cross_sections = photoionization_cross_sections

    def solve(self, electron_temperature, saha_factor, approximation="seaton"):
        """Solve the collisional ionization and recombination rates.

        Parameters
        ----------
        electron_temperature : u.Quantity
            Electron temperatures per cell
        saha_factor : pandas.DataFrame, dtype float
            The Saha factor for each cell. Indexed by atom number, ion number, level number.
        approximation : str, optional
            The rate approximation to use, by default "seaton"

        Returns
        -------
        pd.DataFrame
            Collisional ionization rates
        pd.DataFrame
            Collisional recombination rates

        Raises
        ------
        ValueError
            If an unsupported approximation is requested.
        """
        if approximation == "seaton":
            strength_solver = CollisionalIonizationSeaton(
                self.photoionization_cross_sections
            )
        else:
            raise ValueError(f"approximation {approximation} not supported")

        collision_ionization_rates = strength_solver.solve(electron_temperature)

        # Inverse of the ionization rate for equilibrium
        collision_recombination_rates = collision_ionization_rates.multiply(
            saha_factor.loc[collision_ionization_rates.index]
        )

        return collision_ionization_rates, collision_recombination_rates
