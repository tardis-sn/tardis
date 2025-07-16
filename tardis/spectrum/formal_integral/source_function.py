

from tardis import constants as const

from dataclasses import dataclass
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
from astropy import units as u


class SourceFunctionSolver:

    def __init__(self, line_interaction_type, atomic_data):

        self.line_interaction_type = line_interaction_type
        self.atomic_data = atomic_data


    def solve(self, sim_state, opacity_state, transport_state, levels):
        """
        Calculates the source function using the line absorption rate estimator `Edotlu_estimator`

        Formally it calculates the expression ( 1 - exp(-tau_ul) ) S_ul but this product is what we need later,
        so there is no need to factor out the source function explicitly.
        Parameters
        ----------
        sim_state : tardis.model.SimulationState
        opacity_state : tardis.transport.montecarlo.OpacityState
        transport_state : tardis.transport.montecarlo.TransportState
        atomic_data : tardis.atomic.AtomicData
        levels : pandas.DataFrame
            DataFrame containing the levels of the atomic data

        Returns
        -------
        Numpy array containing ( 1 - exp(-tau_ul) ) S_ul ordered by wavelength of the transition u -> l

        att_S_ul : np.ndarray
            Attenuated source function
        Jredlu : np.ndarray
        Jbluelu : np.ndarray
        e_dot_u : pd.DataFrame
        """

        # Parse states for required values
        v_inner_boundary_index = sim_state.geometry.v_inner_boundary_index
        v_outer_boundary_index = sim_state.geometry.v_outer_boundary_index
        no_of_shells = sim_state.no_of_shells
        dilution_factor = sim_state.dilution_factor
        time_explosion = sim_state.time_explosion
        volume = sim_state.volume

        tau_sobolev = opacity_state.tau_sobolev
        transition_probabilities = opacity_state.transition_probabilities

        j_blue_estimator = transport_state.radfield_mc_estimators.j_blue_estimator
        Edotlu_estimator = transport_state.radfield_mc_estimators.Edotlu_estimator
        time_of_simulation = transport_state.packet_collection.time_of_simulation * u.s


        # slice for the active shells
        local_slice = slice(v_inner_boundary_index,v_outer_boundary_index)

        transition_probabilities = transition_probabilities[:, local_slice]
        tau_sobolevs = tau_sobolev[:, local_slice]

        macro_ref = self.atomic_data.macro_atom_references
        macro_data = self.atomic_data.macro_atom_data

        # levels = self.atomic_data.levels

        no_lvls = len(levels)
        no_shells = len(dilution_factor)

        # Calculate e_dot_u
        upper_level_index = self.atomic_data.lines.index.droplevel(
            "level_number_lower"
        )
        e_dot_u = self.calculate_edotu(time_of_simulation, volume, tau_sobolevs, Edotlu_estimator, 
                        macro_data, macro_ref,
                        transition_probabilities, upper_level_index, no_of_shells, no_lvls,
                        line_interaction_type=self.line_interaction_type
                        )


        # Calculate att_S_ul
        transition_type = self.atomic_data.macro_atom_data.transition_type
        transitions = self.atomic_data.macro_atom_data[transition_type == -1].copy()
        transitions_index = transitions.set_index(
            ["atomic_number", "ion_number", "source_level_number"]
        ).index.copy()
        transition_line_id = transitions.transition_line_id.values
        lines = self.atomic_data.lines.set_index('line_id') 
        lines_idx = lines.index.values 
        
        att_S_ul = self.calculate_att_S_ul(lines, transition_probabilities, 
                            no_of_shells, transition_line_id, lines_idx, 
                            transitions_index, transition_type,
                            e_dot_u, time_explosion)

        # Calculate Jredlu and Jbluelu
        Jbluelu = self.calculate_Jbluelu(time_explosion, time_of_simulation, volume, j_blue_estimator)
        Jredlu = self.calculate_Jredlu(Jbluelu, tau_sobolevs, att_S_ul)

        return SourceFunctionState(att_S_ul, Jredlu, Jbluelu, e_dot_u)


    def calculate_edotu(self, time_of_simulation, volume, tau_sobolevs, Edotlu_estimator, 
                        macro_data, macro_ref,
                        transition_probabilities, upper_level_idx, no_of_shells, no_lvls,
                        line_interaction_type="macroatom"):
        """
        Calculate e_dot_u, the rate energy density is add to the upper level of transitions excited to it
        
        Parameters
        ----------
        time_of_simulation: float
            Time duration of the simulation
        volume: astropy.units.Quantity
        tau_sobolevs: np.ndarray
            Sobolev optical depths 
        Edotlu_estimator: np.ndarray
            The line estimator for the rate of energy absorption of a transition from lower to upper level
        macro_data: pd.DataFrame
            DataFrame containing macro atom data
        macro_ref: pd.DataFrame
            DataFrame containing macro atom references, see http://tardis.readthedocs.io/en/latest/physics/plasma/macroatom.html
        transition_probabilities: np.ndarray
        upper_level_idx: pd.Index
            Index of the upper levels in the atomic data
        no_of_shells: int
            Number of shells in the simulation
        no_lvls: int
            Number of levels in the atomic data
        """

        Edotlu_norm_factor = 1 / (time_of_simulation * volume)
        exptau = 1 - np.exp(-tau_sobolevs)
        Edotlu = Edotlu_norm_factor * exptau * Edotlu_estimator


        columns = range(no_of_shells)
        e_dot_lu = pd.DataFrame(
            Edotlu.value, index=upper_level_idx, columns=columns
        )
        e_dot_u = e_dot_lu.groupby(level=[0, 1, 2]).sum()
        e_dot_u_src_idx = macro_ref.loc[e_dot_u.index].references_idx.values

        e_dot_u.index.names = [
            "atomic_number",
            "ion_number",
            "source_level_number",
        ]

        if line_interaction_type == "macroatom":
            internal_jump_mask = (macro_data.transition_type >= 0).values
            ma_int_data = macro_data[internal_jump_mask]
            internal = transition_probabilities[internal_jump_mask]

            source_level_idx = ma_int_data.source_level_idx.values
            destination_level_idx = ma_int_data.destination_level_idx.values

            C_frame = pd.DataFrame(columns=columns, index=macro_ref.index)
            q_indices = (source_level_idx, destination_level_idx)
            for shell in columns:
                Q = sp.coo_matrix(
                    (internal[:, shell], q_indices), shape=(no_lvls, no_lvls)
                )
                inv_N = sp.identity(no_lvls) - Q
                e_dot_u_vec = np.zeros(no_lvls)
                e_dot_u_vec[e_dot_u_src_idx] = e_dot_u[shell].values
                C_frame[shell] = linalg.spsolve(inv_N.T, e_dot_u_vec)
            
            e_dot_u = C_frame.loc[e_dot_u.index]

        return e_dot_u


    def calculate_att_S_ul(self, lines, transition_probabilities, 
                            no_of_shells, transition_line_id, line_idx, 
                            transitions_index, transition_type,
                            e_dot_u, time_explosion):
        """
        Calculate the attenuated source function

        Parameters
        ----------
        lines: pd.DataFrame
            atomic line data
        transition_probabilities: np.ndarray
        no_of_shells: int
            Number of shells in the simulation
        transition_line_id: np.ndarray
            Line ids for the transitions
        line_idx: np.ndarray
            Indices of the lines in the atomic data
        transitions_index: pd.Index
            Index of the transitions in the macro atom data
        transition_type: np.ndarray
            transition types, see https://tardis-sn.github.io/tardis/physics/setup/plasma/macroatom.html#macroatom for flag definitions
        e_dot_u: pd.DataFrame
            the rate energy density is add to the upper level of transitions excited to it
        time_explosion: float
            geometrical explosion time
        """
        
        q_ul = pd.DataFrame(
            transition_probabilities[
                (transition_type == -1).values
            ],
            index=transitions_index
        )
        wave = lines.wavelength_cm.loc[transition_line_id].values.reshape(-1, 1)
        att_S_ul = wave * (q_ul * e_dot_u) * time_explosion / (4 * np.pi)


        columns = range(no_of_shells)

        result = pd.DataFrame(
            att_S_ul.values,
            index=transition_line_id,
            columns=columns,
        )
        att_S_ul = result.loc[line_idx].values

        return att_S_ul


    def calculate_Jbluelu(self, time_explosion, time_of_simulation, volume, j_blue_estimator):
        """
        Calculates Jbluelu, the normalized J estimator from the blue end of the line from lower to upper level

        Parameters
        ----------
        time_explosion: float
            Time duration of the explosion in seconds
        time_of_simulation: float
            Time duration of the simulation
        volume: astropy.units.Quantity
        j_blue_estimator: np.ndarray
            the line estimator 
        """

        Jbluelu_norm_factor = (
            (
                const.c.cgs
                * time_explosion
                / (
                    4
                    * np.pi
                    * time_of_simulation
                    * volume
                )
            ).to("1/(cm^2 s)").value
        )
        
        # Jbluelu should already by in the correct order, i.e. by wavelength of
        # the transition l->u
        Jbluelu = j_blue_estimator * Jbluelu_norm_factor
        return Jbluelu


    def calculate_Jredlu(self, Jbluelu, tau_sobolevs, att_S_ul):
        """
        Calculates Jredlu, J estimator from the red end of the line from lower to upper level
        """
        
        return Jbluelu * np.exp(-tau_sobolevs) + att_S_ul   


@dataclass
class SourceFunctionState:
    """
    Data class to hold the computed source function values
    """

    att_S_ul: np.ndarray
    Jred_lu: np.ndarray
    Jblue_lu: np.ndarray
    e_dot_u: pd.DataFrame

    def __init__(self, att_S_ul, Jred_lu, Jblue_lu, e_dot_u):
        self.att_S_ul = att_S_ul
        self.Jred_lu = Jred_lu
        self.Jblue_lu = Jblue_lu
        self.e_dot_u = e_dot_u