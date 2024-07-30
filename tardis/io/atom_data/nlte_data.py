import logging

import numpy as np
from scipy import interpolate

logger = logging.getLogger(__name__)


class NLTEData:
    def __init__(self, atom_data, nlte_species):
        self.atom_data = atom_data
        self.lines = atom_data.lines.reset_index()
        self.nlte_species = nlte_species

        if nlte_species:
            logger.info("Preparing the NLTE data")
            self._init_indices()
            if atom_data.collision_data is not None:
                self._create_collision_coefficient_matrix()

    def _init_indices(self):
        self.lines_idx = {}
        self.lines_level_number_lower = {}
        self.lines_level_number_upper = {}
        self.A_uls = {}
        self.B_uls = {}
        self.B_lus = {}

        for species in self.nlte_species:
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

            self.A_uls[species] = self.atom_data.lines.A_ul.values[lines_idx]
            self.B_uls[species] = self.atom_data.lines.B_ul.values[lines_idx]
            self.B_lus[species] = self.atom_data.lines.B_lu.values[lines_idx]

    def _create_collision_coefficient_matrix(self):
        self.C_ul_interpolator = {}
        self.delta_E_matrices = {}
        self.g_ratio_matrices = {}
        collision_group = self.atom_data.collision_data.groupby(
            level=["atomic_number", "ion_number"]
        )
        for species in self.nlte_species:
            no_of_levels = self.atom_data.levels.loc[species].energy.count()
            C_ul_matrix = np.zeros(
                (
                    no_of_levels,
                    no_of_levels,
                    len(self.atom_data.collision_data_temperatures),
                )
            )
            delta_E_matrix = np.zeros((no_of_levels, no_of_levels))
            g_ratio_matrix = np.zeros((no_of_levels, no_of_levels))

            for (
                (
                    atomic_number,
                    ion_number,
                    level_number_lower,
                    level_number_upper,
                ),
                line,
            ) in collision_group.get_group(species).iterrows():
                # line.columns : delta_e, g_ratio, temperatures ...
                C_ul_matrix[
                    level_number_lower, level_number_upper, :
                ] = line.values[2:]
                delta_E_matrix[level_number_lower, level_number_upper] = line[
                    "delta_e"
                ]
                # TODO TARDISATOMIC fix change the g_ratio to be the otherway round - I flip them now here.
                g_ratio_matrix[level_number_lower, level_number_upper] = (
                    1 / line["g_ratio"]
                )
            self.C_ul_interpolator[species] = interpolate.interp1d(
                self.atom_data.collision_data_temperatures, C_ul_matrix
            )
            self.delta_E_matrices[species] = delta_E_matrix

            self.g_ratio_matrices[species] = g_ratio_matrix

    def get_collision_matrix(self, species, t_electrons):
        """
        Creat collision matrix by interpolating the C_ul values for
        the desired temperatures.
        """
        c_ul_matrix = self.C_ul_interpolator[species](t_electrons)
        no_of_levels = c_ul_matrix.shape[0]
        c_ul_matrix[np.isnan(c_ul_matrix)] = 0.0

        # TODO in tardisatomic the g_ratio is the other way round - here I'll flip it in prepare_collision matrix

        c_lu_matrix = (
            c_ul_matrix
            * np.exp(
                -self.delta_E_matrices[species].reshape(
                    (no_of_levels, no_of_levels, 1)
                )
                / t_electrons.reshape((1, 1, t_electrons.shape[0]))
            )
            * self.g_ratio_matrices[species].reshape(
                (no_of_levels, no_of_levels, 1)
            )
        )
        return c_ul_matrix + c_lu_matrix.transpose(1, 0, 2)


class NLTEMatrixFactory:
    def __init__(self, nlte_species):
        pass

    def generate_indices(self):
        pass
