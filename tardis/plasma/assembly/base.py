import pandas as pd
from tardis.util.base import species_string_to_tuple


class PlasmaSolverFactory:

    nlte_species: list
    continuum_interaction_species: pd.MultiIndex

    def __init__(self, config) -> None:
        self.set_nlte_species_from_string(config.plasma.nlte.species)
        self.set_continuum_interaction_species_from_string(
            config.plasma.continuum_interaction.species
        )

    # Convert the continuum interaction species list to a proper format.

    def set_continuum_interaction_species_from_string(
        self, continuum_interaction_species
    ):
        continuum_interaction_species = [
            species_string_to_tuple(species)
            for species in continuum_interaction_species
        ]
        self.continuum_interaction_species = pd.MultiIndex.from_tuples(
            continuum_interaction_species, names=["atomic_number", "ion_number"]
        )

    def set_nlte_species_from_string(self, nlte_species):
        """
        Sets the non-LTE species from a string representation.

        Parameters
        ----------
        nlte_species : str
            A string representation of the non-LTE species.

        Returns
        -------
        None
            This method does not return anything.
        """
        self.nlte_species = [
            species_string_to_tuple(species) for species in nlte_species
        ]
