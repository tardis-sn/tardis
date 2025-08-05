import logging
from dataclasses import dataclass
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
from astropy import units as u

from tardis import constants as const

logger = logging.getLogger(__name__)


class SourceFunctionSolver:
    def __init__(self, line_interaction_type, atomic_data):
        """
        Configures the source function solver

        Parameters
        ----------
        line_interaction_type : str
            The type of line interaction (e.g. "downbranch", "macroatom").
        atomic_data : tardis.atomic.AtomicData
            The atomic data for the simulation.
        """

        self.line_interaction_type = line_interaction_type
        self.atomic_data = atomic_data

    def solve(self, sim_state, opacity_state, transport_state, levels):
        """
        Solves for att_S_ul, Jred_lu, Jblue_lu, and e_dot_u.
        
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
        SourceFunctionState with attributes
            att_S_ul : np.ndarray
                The attenuated source function
            Jred_lu : np.ndarray
                the normalized J estimator from the red end of the line from lower to upper level
            Jblue_lu : np.ndarray
                the normalized J estimator from the blue end of the line from lower to upper level 
            e_dot_u : pd.DataFrame
                The rate energy density is added to the upper level of transitions excited to it
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

        j_blue_estimator = (
            transport_state.radfield_mc_estimators.j_blue_estimator
        )
        e_dot_lu_estimator = (
            transport_state.radfield_mc_estimators.Edotlu_estimator
        )
        time_of_simulation = (
            transport_state.packet_collection.time_of_simulation * u.s
        )

        # slice for the active shells
        local_slice = slice(v_inner_boundary_index, v_outer_boundary_index)

        transition_probabilities = transition_probabilities[:, local_slice]
        tau_sobolevs = tau_sobolev[:, local_slice]

        if self.line_interaction_type == 'macroatom':
            macro_ref = self.atomic_data.macro_atom_references
            macro_data = self.atomic_data.macro_atom_data
        else:
            macro_ref = None
            macro_data = None

        no_lvls = len(levels)
        no_shells = len(dilution_factor)

        # Calculate e_dot_u
        upper_level_index = self.atomic_data.lines.index.droplevel(
            "level_number_lower"
        )
        e_dot_u = self.calculate_e_dot_u(
            time_of_simulation,
            volume,
            tau_sobolevs,
            e_dot_lu_estimator,
            transition_probabilities,
            upper_level_index,
            no_shells,
            no_lvls,
            line_interaction_type=self.line_interaction_type,
            macro_data=macro_data,
            macro_ref=macro_ref
        )

        # Calculate att_S_ul
        transition_type = self.atomic_data.macro_atom_data.transition_type
        transitions = self.atomic_data.macro_atom_data[
            transition_type == -1
        ].copy()
        transitions_index = transitions.set_index(
            ["atomic_number", "ion_number", "source_level_number"]
        ).index.copy()
        transition_line_id = transitions.transition_line_id.values
        lines = self.atomic_data.lines.set_index("line_id")
        lines_idx = lines.index.values

        att_S_ul = self.calculate_att_S_ul(
            lines,
            transition_probabilities,
            no_of_shells,
            transition_line_id,
            lines_idx,
            transitions_index,
            transition_type,
            e_dot_u,
            time_explosion,
        )

        # Calculate Jred_lu and Jblue_lu
        Jblue_lu = self.calculate_Jblue_lu(
            time_explosion, time_of_simulation, volume, j_blue_estimator
        )
        Jred_lu = self.calculate_Jred_lu(Jblue_lu, tau_sobolevs, att_S_ul)

        return SourceFunctionState(att_S_ul, Jred_lu, Jblue_lu, e_dot_u)

    def calculate_e_dot_u(
        self,
        time_of_simulation,
        volume,
        tau_sobolevs,
        e_dot_lu_estimator,
        transition_probabilities,
        upper_level_idx,
        no_of_shells,
        no_lvls,
        line_interaction_type,
        macro_data=None,
        macro_ref=None,
    ):
        """
        Calculate e_dot_u, the rate energy density is added to the upper level of transitions excited to it

        Parameters
        ----------
        time_of_simulation: float
            Time duration of the simulation
        volume: astropy.units.Quantity
        tau_sobolevs: np.ndarray
            Sobolev optical depths
        e_dot_lu_estimator: np.ndarray
            The line estimator for the rate of energy absorption of a transition from lower to upper level
        transition_probabilities: np.ndarray
        upper_level_idx: pd.Index
            Index of the upper levels in the atomic data
        no_of_shells: int
            Number of shells in the simulation
        no_lvls: int
            Number of levels in the atomic data
        line_interaction_type: str
            Type of line interaction (e.g. "macroatom", "downbranch")
        macro_data: pd.DataFrame, optional
            DataFrame containing macro atom data
        macro_ref: pd.DataFrame, optional
            DataFrame containing macro atom references, see http://tardis.readthedocs.io/en/latest/physics/plasma/macroatom.html
        """
        
        # check if macroatom, that the macro atom data exists
        if line_interaction_type == "macroatom" & ((macro_data is None) or (macro_ref is None)):
            raise ValueError(
                "Macro atom data is required for line interaction type 'macroatom'."
            )
        
        if line_interaction_type != "macroatom" & ((macro_data is not None) or (macro_ref is not None)):
            logger.warning(
                "Macro atom data is provided but line interaction type is not 'macroatom'. It will be ignored."
            )

        e_dot_lu_norm_factor = 1 / (time_of_simulation * volume)
        exptau = 1 - np.exp(-tau_sobolevs)
        e_dot_lu = e_dot_lu_norm_factor * exptau * e_dot_lu_estimator

        columns = range(no_of_shells)
        e_dot_lu = pd.DataFrame(
            e_dot_lu.value, index=upper_level_idx, columns=columns
        )
        e_dot_u = e_dot_lu.groupby(level=[0, 1, 2]).sum()
        e_dot_u_src_idx = macro_ref.loc[e_dot_u.index].references_idx.values

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

        # needed for att_S_ul calculation
        e_dot_u.index.names = [
            "atomic_number",
            "ion_number",
            "source_level_number",
        ]

        return e_dot_u

    def calculate_att_S_ul(
        self,
        lines,
        transition_probabilities,
        no_of_shells,
        transition_line_id,
        line_idx,
        transitions_index,
        transition_type,
        e_dot_u,
        time_explosion,
    ):
        """
        Calculates the source function using the line absorption rate estimator `e_dot_lu_estimator`

        Formally it calculates the expression ( 1 - exp(-tau_ul) ) S_ul but this product is what we need later,
        so there is no need to factor out the source function explicitly.

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
            transition_probabilities[(transition_type == -1).values],
            index=transitions_index,
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

    def calculate_Jblue_lu(
        self, time_explosion, time_of_simulation, volume, j_blue_estimator
    ):
        """
        Calculates Jblue_lu, the normalized J estimator from the blue end of the line from lower to upper level

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

        Jblue_lu_norm_factor = (
            (
                const.c.cgs
                * time_explosion
                / (4 * np.pi * time_of_simulation * volume)
            )
            .to("1/(cm^2 s)")
            .value
        )

        # Jblue_lu should already by in the correct order, i.e. by wavelength of
        # the transition l->u
        Jblue_lu = j_blue_estimator * Jblue_lu_norm_factor
        return Jblue_lu

    def calculate_Jred_lu(self, Jblue_lu, tau_sobolevs, att_S_ul):
        """
        Calculates Jred_lu, J estimator from the red end of the line from lower to upper level
        """

        return Jblue_lu * np.exp(-tau_sobolevs) + att_S_ul


@dataclass
class SourceFunctionState:
    """
    Data class to hold the computed source function values

    Attributes:
    ----------
    
    att_S_ul : np.ndarray
        The attenuated source function
    Jred_lu : np.ndarray
        the normalized J estimator from the red end of the line from lower to upper level
    Jblue_lu : np.ndarray
        the normalized J estimator from the blue end of the line from lower to upper level 
    e_dot_u : pd.DataFrame
        The rate energy density is added to the upper level of transitions excited to it
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
