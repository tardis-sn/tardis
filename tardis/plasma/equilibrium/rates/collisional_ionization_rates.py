from tardis.plasma.equilibrium.rates.collisional_ionization_strengths import (
    CollisionalIonizationSeaton,
)


class CollisionalIonizationRateSolver:
    """Solver for collisional ionization and recombination rates."""

    def __init__(self, photoionization_cross_sections):
        """Initialize the collisional ionization rate solver.

        Parameters
        ----------
        photoionization_cross_sections : pd.DataFrame
            Photoionization cross sections.
        """
        self.photoionization_cross_sections = photoionization_cross_sections

    @staticmethod
    def __reindex_ionization_rate_dataframe(
        rate_dataframe, recombination=False
    ):
        """Index the ionization rate dataframe to include source and destination
        ion numbers and level numbers.

        Parameters
        ----------
        rate_dataframe : pd.DataFrame
            Dataframe of ionization rates
        recombination : bool, optional
            If true, reverse the direction of source to destination
            to handle recombination, by default False

        Returns
        -------
        pd.DataFrame
            Dataframe with additional columns for source and destination
            ion numbers and level numbers, indexed by atomic number, ion number,
            source level number, destination ion number, and destination
            level number.
        """
        rate_dataframe.index.names = [
            "atomic_number",
            "ion_number",
            "level_number_source",
        ]

        rate_dataframe = rate_dataframe.reset_index()

        if recombination:
            rate_dataframe["ion_number_destination"] = rate_dataframe[
                "ion_number"
            ]
            rate_dataframe["ion_number_source"] = (
                rate_dataframe["ion_number"] + 1
            )
        else:
            rate_dataframe["ion_number_source"] = rate_dataframe["ion_number"]
            rate_dataframe["ion_number_destination"] = (
                rate_dataframe["ion_number"] + 1
            )

        # ionized electrons are assumed to leave the ion in the ground state for now
        rate_dataframe["level_number_destination"] = 0

        not_fully_ionized_mask = (
            rate_dataframe["atomic_number"] != rate_dataframe["ion_number"]
        )

        rate_dataframe = rate_dataframe[not_fully_ionized_mask]

        rate_dataframe = rate_dataframe.set_index(
            [
                "atomic_number",
                "ion_number",
                "ion_number_source",
                "ion_number_destination",
                "level_number_source",
                "level_number_destination",
            ]
        )

        return rate_dataframe

    def solve(self, electron_distribution, saha_factor, approximation="seaton"):
        """Solve the collisional ionization and recombination rates.

        Parameters
        ----------
        electron_distribution :
            Electron distribution per cell
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

        collision_ionization_rates = strength_solver.solve(
            electron_distribution.temperature
        )

        # Inverse of the ionization rate for equilibrium
        collision_recombination_rates = collision_ionization_rates.multiply(
            saha_factor
        )

        collision_ionization_rates = (
            self.__reindex_ionization_rate_dataframe(
                collision_ionization_rates, recombination=False
            )
            * electron_distribution.number_density
        )

        collision_recombination_rates = (
            self.__reindex_ionization_rate_dataframe(
                collision_recombination_rates, recombination=True
            )
        ) * electron_distribution.number_density**2

        return collision_ionization_rates, collision_recombination_rates
