

from tardis import constants as const
from tardis.spectrum.formal_integral.base import interpolate_integrator_quantities

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as linalg

from astropy import units as u


class SourceFunctionSolver:

    def __init__(self, line_interaction_type):

        # self.configuration = configuration
        self.line_interaction_type = line_interaction_type


    def solve(self, sim_state, opacity_state, transport_state, atomic_data, levels):
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

        # TODO: check if the opacity state only lives in the transport state
            # and if so, remove the opacity_state parameter
            # and use transport_state.opacity_state instead
        tau_sobolev = opacity_state.tau_sobolev
        transition_probabilities = opacity_state.transition_probabilities

        j_blue_estimator = transport_state.radfield_mc_estimators.j_blue_estimator
        Edotlu_estimator = transport_state.radfield_mc_estimators.Edotlu_estimator
        time_of_simulation = transport_state.packet_collection.time_of_simulation * u.s


        # slice for the active shells
        local_slice = slice(v_inner_boundary_index,v_outer_boundary_index)

        transition_probabilities = transition_probabilities[:, local_slice]
        tau_sobolevs = tau_sobolev[:, local_slice]

        columns = range(no_of_shells)

        macro_ref = atomic_data.macro_atom_references
        macro_data = atomic_data.macro_atom_data

        no_lvls = len(levels)
        no_shells = len(dilution_factor) # TODO: is dilution_factor len=no_of_shells 

        # Calculate e_dot_u
        upper_level_index = atomic_data.lines.index.droplevel(
            "level_number_lower"
        )
        e_dot_u = calculate_edotu(time_of_simulation, volume, tau_sobolevs, Edotlu_estimator, 
                        macro_data, macro_ref, self.line_interaction_type,
                        transition_probabilities, upper_level_index, columns, no_lvls)


        # Calculate att_S_ul
        transitions_index = transitions.set_index(
            ["atomic_number", "ion_number", "source_level_number"]
        ).index.copy()
        transition_type = atomic_data.macro_atom_data.transition_type
        transitions = atomic_data.macro_atom_data[transition_type == -1].copy()
        transition_line_id = transitions.transition_line_id.values
        lines = atomic_data.lines.set_index('line_id') # TODO: investigate why this is like this
        lines_idx = lines.index.values 
        
        att_S_ul = calculate_att_S_ul(lines, transition_probabilities, 
                            no_of_shells, transition_line_id, lines_idx, 
                            transitions_index, transition_type,
                            e_dot_u, time_explosion)

        # Calculate Jredlu and Jbluelu
        Jbluelu = calculate_Jbluelu(time_explosion, time_of_simulation, volume, j_blue_estimator)
        Jredlu = calculate_Jredlu(Jbluelu, tau_sobolevs, att_S_ul)

        return att_S_ul, Jredlu, Jbluelu, e_dot_u
        
    
# transport
    # time of sim, Edotlu_estimator, line_interaction_type
# sim state
    # volume
# upper_level_index = atomic_data.lines.index.droplevel("level_number_lower")
def calculate_edotu(time_of_simulation, volume, tau_sobolevs, Edotlu_estimator, 
                    macro_data, macro_ref, line_interaction_type,
                    transition_probabilities, upper_level_idx, columns, no_lvls):

    Edotlu_norm_factor = 1 / (time_of_simulation * volume)
    exptau = 1 - np.exp(-tau_sobolevs)
    Edotlu = Edotlu_norm_factor * exptau * Edotlu_estimator

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

    # if transport.line_interaction_type == "macroatom":
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
    # if transport.line_interaction_type == "macroatom":
        e_dot_u = C_frame.loc[e_dot_u.index]

    return e_dot_u


# transition_line_id: 
    # transitions = macro_atom_data[macro_atom_data.transition_type == -1].copy()
    # transition_line_id = transitions.transition_line_id.values
# line_idx:
    # line_idx = 
    # line_idx = lines.index.values
def calculate_att_S_ul(lines, transition_probabilities, 
                        no_of_shells, transition_line_id, line_idx, 
                        transitions_index, transition_type,
                        e_dot_u, time_explosion):
    
    # rewrite as : q_ul = dataframe(index = trans_idx)
    q_ul = pd.DataFrame(
        transition_probabilities[
            (transition_type == -1).values
        ], 
        index = transitions_index
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

# time_explosion and volume from sim state
# time of sim from montecarlo transport state which is transport state
# j_blue_estimator from transport state
def calculate_Jbluelu(time_explosion, time_of_simulation, volume, j_blue_estimator):

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


def calculate_Jredlu(Jbluelu, tau_sobolevs, att_S_ul):
    
    return Jbluelu * np.exp(-tau_sobolevs) + att_S_ul   