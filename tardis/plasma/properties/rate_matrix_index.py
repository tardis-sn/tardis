import pandas as pd
from tardis.plasma.properties.base import ProcessingPlasmaProperty

__all__ = [
    "NLTEIndexHelper",
]


class NLTEIndexHelper(ProcessingPlasmaProperty):
    outputs = ("rate_matrix_index",)

    def __init__(self, plasma_parent, nlte_ionization_species=0):
        super().__init__(plasma_parent)
        self.nlte_ionization_species = nlte_ionization_species

    def calculate(self, levels, nlte_ionization_species):
        nlte_excitation_species = []  # not yet implemented
        rate_matrix_index = pd.MultiIndex.from_tuples(
            list(
                self.calculate_rate_matrix_index(
                    levels,
                    self.nlte_ionization_species,
                    nlte_excitation_species,
                )
            ),
            names=levels.names,
        ).drop_duplicates()
        return rate_matrix_index

    def calculate_rate_matrix_index(self, levels, nlte_excitation_species=[]):
        for level in levels:
            if level[:2] in self.nlte_ionization_species:
                yield (*level[:2], "nlte_ion")
            elif (level[:2] not in self.nlte_ionization_species) and (
                level[:2] not in nlte_excitation_species
            ):
                yield (*level[:2], "lte_ion")
            else:
                yield level
        yield ("n_e", "n_e", "n_e")
