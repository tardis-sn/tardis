import pandas as pd


class ContinuumInputData(object):
    """
    The common input data object for all continuum calculations.
    """
    def __init__(self, atom_data, plasma_array, ws, radiative_transition_probabilities):
        # Plasma quantities
        self.electron_densities = plasma_array.electron_densities.values
        self.t_electrons = plasma_array.t_electrons
        self.t_rads = plasma_array.t_rad
        self.link_t_rad_t_electron = plasma_array.link_t_rad_t_electron
        self.ion_number_density = plasma_array.ion_number_density
        # TODO: Replace level population with LTE level population
        self.lte_level_pop = plasma_array.level_number_density
        self.level_pop = plasma_array.level_number_density

        # Radiation field
        self.ws = ws

        # Atom data
        self.lines = atom_data.lines
        self.levels = atom_data.levels
        self.ionization_energies = atom_data.ionization_data
        self.photoionization_data = atom_data.continuum_data.photoionization_data

        #
        self.macro_atom_references = atom_data.macro_atom_references
        self.macro_atom_references_by_idx = atom_data.macro_atom_references.reset_index().set_index('references_idx')
        self.macro_atom_data = atom_data.macro_atom_data
        self.macro_atom_continuum_data = atom_data.continuum_data.macro_atom_data
        self.radiative_transition_probabilities = radiative_transition_probabilities
        self.radiative_transition_probabilities_prep = self._prepare_radiative_probabilities(
            radiative_transition_probabilities)

        #
        self.continuum_references = atom_data.continuum_data.continuum_references
        self.continuum_data = atom_data.continuum_data.continuum_data

        #
        self.ion_charges = self._get_ion_charges()

        # Computed quantities
        self.nu_i = self._get_nu_i()

    def _get_nu_i(self):
        return self.photoionization_data.groupby(level=[0, 1, 2]).first().nu.values

    def _get_ion_charges(self):
        return self.ion_number_density.index.get_level_values(1).values

    def _prepare_radiative_probabilities(self, radiative_prob):
        source_level_idx = self._get_source_level_idx()
        destination_level_idx = self.macro_atom_data.destination_level_idx.values
        new_index = pd.MultiIndex.from_arrays([source_level_idx, destination_level_idx],
                                              names=['source_level_idx', 'destination_level_idx'])
        radiative_prob_prep = radiative_prob.set_index(new_index)
        return radiative_prob_prep

    def _get_source_level_idx(self):
        macro_atom_data = self.macro_atom_data
        source_level_index = pd.MultiIndex.from_arrays([macro_atom_data['atomic_number'], macro_atom_data['ion_number'],
                                                        macro_atom_data['source_level_number']])
        return self._get_level_idx(source_level_index)

    def _get_level_idx(self, multi_index):
        return self.macro_atom_references.loc[multi_index, 'references_idx'].values