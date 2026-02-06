import numpy as np
import pandas as pd


class ContinuumInputData:
    """
    The common input data object for all continuum calculations.
    """

    def __init__(
        self,
        atom_data,
        plasma_array,
        ws,
        radiative_transition_probabilities,
        estimators,
    ):
        self.selected_continuum_species = plasma_array.nlte_species
        # Plasma quantities
        self.electron_densities = plasma_array.electron_densities.values
        self.t_electrons = plasma_array.t_electrons
        self.t_rads = plasma_array.t_rad
        self.link_t_rad_t_electron = plasma_array.link_t_rad_t_electron
        self.ion_number_density = plasma_array.ion_number_density
        self.lte_ion_number_density = plasma_array.lte_ion_number_density
        self.lte_level_pop = plasma_array.lte_level_number_density
        self.level_pop = plasma_array.level_number_density
        self.q_ij = plasma_array.coll_exc_coeff
        self.coll_ion_coeff = plasma_array.coll_ion_coeff
        self.coll_recomb_coeff = plasma_array.coll_recomb_coeff
        self.alpha_sp = plasma_array.alpha_sp
        self.sp_fb_cooling_rates = plasma_array.sp_fb_cooling_rates

        # Radiation field
        self.ws = ws
        self.estimators = estimators

        # Atom data
        self.atom_data = atom_data
        # TODO: Unclear if necessary. Is lines used?
        mask = self._get_continuum_species_mask(
            self.atom_data, self.selected_continuum_species
        )
        self.lines = atom_data.lines[mask]
        # self.lines = atom_data.lines[atom_data.lines['atomic_number'].isin(selected_continuum_species)]

        self.levels = atom_data.levels
        self.ionization_energies = atom_data.ionization_data
        self.photoionization_data = (
            atom_data.continuum_data.photoionization_data
        )

        self.macro_atom_references = atom_data.macro_atom_references
        self.macro_atom_references_by_idx = (
            atom_data.macro_atom_references.reset_index().set_index(
                "references_idx"
            )
        )
        self.macro_atom_data = atom_data.macro_atom_data
        # self.macro_atom_continuum_data = atom_data.continuum_data.macro_atom_data
        self.radiative_transition_probabilities = (
            radiative_transition_probabilities
        )
        self.radiative_transition_probabilities_prep = (
            self._prepare_radiative_probabilities(
                radiative_transition_probabilities
            )
        )

        self.continuum_references = (
            atom_data.continuum_data.continuum_references
        )
        self.continuum_data = atom_data.continuum_data.continuum_data

        self.ion_charges = self._get_ion_charges()

        # Computed quantities
        self.no_lvls_cont_species = self._get_no_lvls_continuum_species()
        self.ref_idx_lvls_cont_species = (
            self._get_references_idx_continuum_species_lvls()
        )
        self.is_diag_macro_level = self._get_diag_levels_mask()
        self.diag_lvls_start_idx = self._get_diag_lvls_start_idx(
            atom_data.macro_atom_references, self.selected_continuum_species
        )
        self.diag_lvls_block_ref = self._get_diag_lvls_block_ref(
            atom_data.macro_atom_references, self.selected_continuum_species
        )
        self.nu_i = self._get_nu_i()

    @staticmethod
    def _get_diag_lvls_start_idx(
        macro_atom_references, selected_continuum_species
    ):
        diag_lvls_start_idx = [
            macro_atom_references.loc[species, "references_idx"].iloc[0]
            for species in selected_continuum_species
        ]
        return np.array(diag_lvls_start_idx)

    @staticmethod
    def _get_diag_lvls_block_ref(
        macro_atom_references, selected_continuum_species
    ):
        diag_lvls_block_ref = [
            len(macro_atom_references.loc[species, "references_idx"])
            for species in selected_continuum_species
        ]
        return np.hstack([[0], np.array(diag_lvls_block_ref).cumsum()])

    @staticmethod
    def _get_continuum_species_mask(atom_data, selected_continuum_species):
        atomic_number = atom_data.lines["atomic_number"].values
        ion_number = atom_data.lines["ion_number"].values
        mask = np.zeros_like(atomic_number, dtype=bool)
        for species in selected_continuum_species:
            species_mask = np.logical_and(
                atomic_number == species[0], ion_number == species[1]
            )
            mask = np.logical_or(mask, species_mask)
        return mask

    def _get_nu_i(self):
        nu_i = self.photoionization_data.groupby(level=[0, 1, 2]).first().nu
        return nu_i

    def _get_ion_charges(self):
        return self.ion_number_density.index.get_level_values(1).values

    def _get_no_lvls_continuum_species(self):
        no_lvls = 0
        for species in self.selected_continuum_species:
            no_lvls += len(self.macro_atom_references.loc[species])
        return no_lvls

    def _get_references_idx_continuum_species_lvls(self):
        species_idx_list = []
        for species in self.selected_continuum_species:
            species_idx_list.append(
                self.macro_atom_references.loc[species, "references_idx"].values
            )
        return np.hstack(species_idx_list)

    # Needed in cmontecarlo to test if activate idx corresponds to diag MA
    def _get_diag_levels_mask(self):
        atoms = self.macro_atom_references.index.get_level_values(0).values
        ions = self.macro_atom_references.index.get_level_values(1).values
        mask = np.zeros(len(self.macro_atom_references), dtype=bool)
        for species in self.selected_continuum_species:
            mask_atom = species[0] == atoms
            mask_ion = species[1] == ions
            mask_species = np.logical_and(mask_atom, mask_ion)
            mask = np.logical_or(mask, mask_species)
        return mask.astype(int)

    def _prepare_radiative_probabilities(self, radiative_prob):
        source_level_idx = self._get_source_level_idx()
        destination_level_idx = (
            self.macro_atom_data.destination_level_idx.values
        )
        new_index = pd.MultiIndex.from_arrays(
            [source_level_idx, destination_level_idx],
            names=["source_level_idx", "destination_level_idx"],
        )
        radiative_prob_prep = radiative_prob.set_index(new_index)
        return radiative_prob_prep

    def _get_source_level_idx(self):
        macro_atom_data = self.macro_atom_data
        source_level_index = pd.MultiIndex.from_arrays(
            [
                macro_atom_data["atomic_number"],
                macro_atom_data["ion_number"],
                macro_atom_data["source_level_number"],
            ]
        )
        return self._get_level_idx(source_level_index)

    def _get_level_idx(self, multi_index):
        return self.macro_atom_references.loc[
            multi_index, "references_idx"
        ].values
