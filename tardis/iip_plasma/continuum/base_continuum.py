from tardis.iip_plasma.continuum.base import (
    InverseProcess,
    PhysicalContinuumProcess,
)
from tardis.iip_plasma.continuum.collisional_processes import *
from tardis.iip_plasma.continuum.cooling import CoolingRates
from tardis.iip_plasma.continuum.exceptions import InvalidContinuumProcessError
from tardis.iip_plasma.continuum.input_data import ContinuumInputData
from tardis.iip_plasma.continuum.probabilities import (
    RecombinationTransitionProbabilities,
    TransitionProbabilities,
)
from tardis.iip_plasma.continuum.radiative_processes import *

default_processes = [
    "radiative_recombination",
    "free_free",
    "collisional_recombination",
    "radiative_excitation",
    "collisional_ionization",
    "collisional_excitation",
    "collisional_deexcitation",
    "radiative_deexcitation",
    "radiative_ionization",
]


class BaseContinuum:
    direct_processes = {
        process.name: process
        for process in PhysicalContinuumProcess.__subclasses__()
    }
    inverse_processes = [
        process.name for process in InverseProcess.__subclasses__()
    ]
    process2inverse_process = {
        process.name_of_inverse_process: process
        for process in InverseProcess.__subclasses__()
    }

    def __init__(
        self,
        atom_data,
        plasma_array,
        ws,
        radiative_transition_probabilities,
        estimators,
        requested_processes=default_processes,
    ):
        self._validate_requested_processes(requested_processes)
        self.input = ContinuumInputData(
            atom_data,
            plasma_array,
            ws,
            radiative_transition_probabilities,
            estimators,
        )
        self._set_physical_processes(requested_processes)
        self._set_inverse_processes()
        self._set_cooling_rates()
        self._set_recombination_transition_probabilities()
        self._set_transition_probabilities()
        self._diagonalize_ma()

    def _set_physical_processes(self, requested_processes):
        for name, process in BaseContinuum.direct_processes.iteritems():
            if name in requested_processes:
                setattr(self, name, process(self.input))

    def _set_inverse_processes(self):
        for (
            name_of_inverse_process,
            process,
        ) in BaseContinuum.process2inverse_process.iteritems():
            if hasattr(self, name_of_inverse_process):
                inverse_process = getattr(self, name_of_inverse_process)
                setattr(
                    self,
                    process.name,
                    process.from_inverse_process(inverse_process),
                )

    def _set_cooling_rates(self):
        cooling_processes = {
            name: process
            for name, process in self.__dict__.iteritems()
            if hasattr(process, "cooling") and process.cooling is True
        }
        self.cooling_rates = CoolingRates(self.input, **cooling_processes)

    def _set_transition_probabilities(self):
        macro_atom_processes = {
            name: process
            for name, process in self.__dict__.iteritems()
            if hasattr(process, "macro_atom_transitions")
            and process.macro_atom_transitions is not None
        }
        self.transition_probabilities = TransitionProbabilities(
            self.input, **macro_atom_processes
        )

    def _set_recombination_transition_probabilities(self):
        recombination_process_list = [
            "radiative_recombination",
            "collisional_recombination",
        ]
        recombination_processes = {}
        for process in recombination_process_list:
            if hasattr(self, process):
                recombination_processes[process] = getattr(self, process)
        self.recombination_transition_probabilities = (
            RecombinationTransitionProbabilities(
                self.input, **recombination_processes
            )
        )

    def _get_ref_indices(self):
        x_ind = np.tile(
            self.input.ref_idx_lvls_cont_species.reshape(-1, 1),
            len(self.input.ref_idx_lvls_cont_species),
        )
        y_ind = x_ind.T
        return x_ind, y_ind

    def _diagonalize_ma(self):
        no_shells = len(self.input.t_electrons)
        # no_lvls = len(self.input.macro_atom_references.loc[(1,0)])
        no_lvls = self.input.no_lvls_cont_species

        no_cont_lvls = len(self.input.selected_continuum_species)
        self.no_cont_lvls = no_cont_lvls
        no_lvls_tot = no_lvls + no_cont_lvls + 1
        x_ind, y_ind = self._get_ref_indices()
        index = zip(x_ind.flatten(), y_ind.flatten())
        ref_idx_2_pos = pd.Series(
            np.arange(no_lvls, dtype=np.int64),
            index=self.input.ref_idx_lvls_cont_species,
        )
        no_k_i = no_cont_lvls + 1

        trans_prob = self.transition_probabilities.dataframe
        cont_jump_mask = trans_prob.transition_type == 2
        internal_mask = np.logical_or(
            trans_prob.transition_type == 1, trans_prob.transition_type == 0
        )
        transient_mask = np.logical_or(cont_jump_mask, internal_mask)
        source_level = trans_prob.index.get_level_values("source_level_idx")

        int_jump_prob = self.transition_probabilities.dataframe.query(
            "transition_type == 0 | transition_type == 1"
        )
        cont_jump_prob = self.transition_probabilities.dataframe.query(
            "transition_type == 2"
        )
        cont_ma_jump_prob = (
            self.recombination_transition_probabilities.dataframe.query(
                "transition_type == 0"
            )
        )
        self.B = np.zeros((no_shells, no_lvls_tot - 1, no_lvls_tot - 1))
        self.expected_steps = np.zeros((no_shells, no_lvls_tot - 1))

        for shell in range(no_shells):
            Q = np.zeros((no_lvls_tot, no_lvls_tot))

            # Setup internal MA jump matrix
            int_jump_matrix = (
                int_jump_prob.loc[index, shell]
                .fillna(0)
                .values.reshape((no_lvls, no_lvls))
            )
            Q[:-no_k_i, :-no_k_i] = int_jump_matrix

            # Setup jumps from MA to continuum
            cj_source_idx = cont_jump_prob.index.get_level_values(0).values
            cj_source_idx = ref_idx_2_pos.loc[cj_source_idx].values
            cj_dest_idx = (
                cont_jump_prob.index.get_level_values(1).values + no_lvls
            )
            cj_index = (cj_source_idx, cj_dest_idx)
            Q[cj_index] = cont_jump_prob[shell].values

            # Setup jumps from continuum to MA
            c_ma_dest_idx = cont_ma_jump_prob.destination_level_idx.values
            c_ma_dest_idx = ref_idx_2_pos.loc[c_ma_dest_idx].values

            c_ma_source_idx_old = np.zeros_like(c_ma_dest_idx) + no_lvls
            c_ma_source_idx = (
                cont_jump_prob.index.get_level_values(1).values + no_lvls
            )
            c_ma_index = (c_ma_source_idx, c_ma_dest_idx)
            Q[c_ma_index] = cont_ma_jump_prob[shell].values

            # Calculate the absorption probabilities
            Qs = Q[:-1, :-1]
            inv_N = np.identity(Qs.shape[0]) - Qs
            N = np.linalg.inv(inv_N)
            expected_steps = N.sum(axis=1)
            R = np.diag(1 - Qs.sum(axis=1))
            B = np.dot(N, R)
            self.B[shell] = B
            self.expected_steps[shell] = expected_steps

        cont_lvls_mask = source_level.isin(self.input.ref_idx_lvls_cont_species)
        internal_drop_mask = np.logical_and(transient_mask, cont_lvls_mask)
        trans_prob = trans_prob[np.logical_not(internal_drop_mask)]
        # TODO: Do not renormalize the entire dataframe
        trans_prob = (
            self.transition_probabilities._normalize_transition_probabilities(
                trans_prob, no_ref_columns=2
            )
        )
        self.transition_probabilities.dataframe = trans_prob
        self.transition_probabilities._set_montecarlo_data()

        recomb_prob = (
            self.recombination_transition_probabilities.dataframe.query(
                "transition_type==-2 | transition_type==-3"
            )
        )
        recomb_prob = (
            self.transition_probabilities._normalize_transition_probabilities(
                recomb_prob, no_ref_columns=3
            )
        )
        # TODO: Is this sufficient to guarantee well defined block references
        recomb_prob.sort_index(inplace=True)
        self.recombination_transition_probabilities.dataframe = recomb_prob
        self.recombination_transition_probabilities._set_montecarlo_data()

    @classmethod
    def _validate_requested_processes(cls, requested_processes):
        for process_name in requested_processes:
            if (
                process_name in cls.direct_processes
                or process_name in cls.inverse_processes
            ):
                continue
            else:
                raise InvalidContinuumProcessError(process_name)
