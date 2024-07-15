import logging

import numpy as np

logger = logging.getLogger(__name__)


class NLTEExcitationData(object):
    """Data needed for NLTE excitation treatment."""

    def __init__(self, atom_data_lines, nlte_excitation_species):
        """Initializes the NLTEExcitationData object.

        Parameters
        ----------
        atom_data_lines : pandas.DataFrame
            Information on lines in the atomic data.
        nlte_excitation_species : list
            List of species treated in NLTE excitation.
        """
        self.atom_data_lines = atom_data_lines
        self.lines = atom_data_lines.reset_index()
        self.nlte_excitation_species = nlte_excitation_species

        if nlte_excitation_species:
            logger.info("Preparing the NLTE data")
            self._init_indices()

    def _init_indices(self):
        """Initializes A_ul and B_ul, B_lu coefficients."""
        self.lines_idx = {}
        self.lines_level_number_lower = {}
        self.lines_level_number_upper = {}
        self.A_uls = {}
        self.B_uls = {}
        self.B_lus = {}
        nlte_excitation_species = self.nlte_excitation_species
        for species in nlte_excitation_species:
            lines_idx = np.where(
                (self.lines.atomic_number == species[0])
                & (self.lines.ion_number == species[1])
            )
            self.lines_idx[species] = lines_idx
            self.lines_level_number_lower[
                species
            ] = self.lines.level_number_lower.values[lines_idx].astype(int)
            self.lines_level_number_upper[
                species
            ] = self.lines.level_number_upper.values[lines_idx].astype(int)

            self.A_uls[species] = self.atom_data_lines.A_ul.values[lines_idx]
            self.B_uls[species] = self.atom_data_lines.B_ul.values[lines_idx]
            self.B_lus[species] = self.atom_data_lines.B_lu.values[lines_idx]
