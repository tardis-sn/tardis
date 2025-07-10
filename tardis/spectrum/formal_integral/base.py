import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import scipy.sparse as sp
import scipy.sparse.linalg as linalg

import warnings

from tardis import constants as const
from tardis.transport.montecarlo.configuration import montecarlo_globals

C_INV = 3.33564e-11
KB_CGS = 1.3806488e-16
H_CGS = 6.62606957e-27

# here will go function that are independent of the type of formal integral
class BoundsError(IndexError):
    pass


class IntegrationError(Exception):
    pass


def check(simulation_state, plasma, transport, raises=True):
    """
    A method that determines if the formal integral can be performed with
    the current configuration settings

    The function returns False if the configuration conflicts with the
    required settings. If raises evaluates to True, then a
    IntegrationError is raised instead
    """

    def raise_or_return(message):
        if raises:
            raise IntegrationError(message)
        else:
            warnings.warn(message)
            return False

    for obj in (simulation_state, plasma, transport):
        if obj is None:
            return raise_or_return(
                "The integrator is missing either model, plasma or "
                "transport. Please make sure these are provided to the "
                "FormalIntegrator."
            )

    if transport.line_interaction_type not in [
        "downbranch",
        "macroatom",
    ]:
        return raise_or_return(
            "The FormalIntegrator currently only works for "
            'line_interaction_type == "downbranch"'
            'and line_interaction_type == "macroatom"'
        )

    if montecarlo_globals.CONTINUUM_PROCESSES_ENABLED:
        return raise_or_return(
            "The FormalIntegrator currently does not work for "
            "continuum interactions."
        )

    return True


def calculate_p_values(R_max, N):
    """
    Calculates the p values of N

    Parameters
    ----------
    R_max : float64
    N : int64

    Returns
    -------
    float64
    """
    return np.arange(N).astype(np.float64) * R_max / (N - 1)


def intensity_black_body(nu, temperature):
    """
    Calculate the blackbody intensity.

    Parameters
    ----------
    nu : float64
        frequency
    temperature : float64
        Temperature

    Returns
    -------
    float64
    """
    if nu == 0:
        return np.nan  # to avoid ZeroDivisionError
    beta_rad = 1 / (KB_CGS * temperature)
    coefficient = 2 * H_CGS * C_INV * C_INV
    return coefficient * nu * nu * nu / (np.exp(H_CGS * nu * beta_rad) - 1)


def make_source_function(simulation_state, opacity_state, transport, plasma, interpolate_shells=0):
    """
    Calculates the source function using the line absorption rate estimator `Edotlu_estimator`

    Formally it calculates the expression ( 1 - exp(-tau_ul) ) S_ul but this product is what we need later,
    so there is no need to factor out the source function explicitly.
    Parameters
    ----------
    model : tardis.model.SimulationState

    Returns
    -------
    Numpy array containing ( 1 - exp(-tau_ul) ) S_ul ordered by wavelength of the transition u -> l
    """

    montecarlo_transport_state = transport.transport_state
    atomic_data = plasma.atomic_data
    levels_index = plasma.levels

    # slice for the active shells
    local_slice = slice(
        simulation_state.geometry.v_inner_boundary_index,
        simulation_state.geometry.v_outer_boundary_index,
    )

    transition_probabilities = opacity_state.transition_probabilities[
        :, local_slice
    ]
    tau_sobolevs = opacity_state.tau_sobolev[:, local_slice]

    columns = range(simulation_state.no_of_shells)

    # macro_ref = atomic_data.macro_atom_references
    macro_ref = atomic_data.macro_atom_references
    # macro_data = atomic_data.macro_atom_data
    macro_data = plasma.atomic_data.macro_atom_data

    no_lvls = len(levels_index)
    no_shells = len(simulation_state.dilution_factor)

    if transport.line_interaction_type == "macroatom":
        internal_jump_mask = (macro_data.transition_type >= 0).values
        ma_int_data = macro_data[internal_jump_mask]
        internal = transition_probabilities[internal_jump_mask]

        source_level_idx = ma_int_data.source_level_idx.values
        destination_level_idx = ma_int_data.destination_level_idx.values

    Edotlu_norm_factor = 1 / (
        montecarlo_transport_state.packet_collection.time_of_simulation
        * simulation_state.volume
    )
    exptau = 1 - np.exp(-tau_sobolevs)
    Edotlu = (
        Edotlu_norm_factor
        * exptau
        * montecarlo_transport_state.radfield_mc_estimators.Edotlu_estimator
    )

    # The following may be achieved by calling the appropriate plasma
    # functions
    Jbluelu_norm_factor = (
        (
            const.c.cgs
            * simulation_state.time_explosion
            / (
                4
                * np.pi
                * montecarlo_transport_state.time_of_simulation
                * simulation_state.volume
            )
        )
        .to("1/(cm^2 s)")
        .value
    )
    # Jbluelu should already by in the correct order, i.e. by wavelength of
    # the transition l->u
    Jbluelu = (
        transport.transport_state.radfield_mc_estimators.j_blue_estimator
        * Jbluelu_norm_factor
    )

    upper_level_index = atomic_data.lines.index.droplevel(
        "level_number_lower"
    )
    e_dot_lu = pd.DataFrame(
        Edotlu.value, index=upper_level_index, columns=columns
    )
    e_dot_u = e_dot_lu.groupby(level=[0, 1, 2]).sum()
    e_dot_u_src_idx = macro_ref.loc[e_dot_u.index].references_idx.values

    if transport.line_interaction_type == "macroatom":
        C_frame = pd.DataFrame(columns=columns, index=macro_ref.index)
        q_indices = (source_level_idx, destination_level_idx)
        for shell in range(no_shells):
            Q = sp.coo_matrix(
                (internal[:, shell], q_indices), shape=(no_lvls, no_lvls)
            )
            inv_N = sp.identity(no_lvls) - Q
            e_dot_u_vec = np.zeros(no_lvls)
            e_dot_u_vec[e_dot_u_src_idx] = e_dot_u[shell].values
            C_frame[shell] = linalg.spsolve(inv_N.T, e_dot_u_vec)

    e_dot_u.index.names = [
        "atomic_number",
        "ion_number",
        "source_level_number",
    ]  # To make the q_ul e_dot_u product work, could be cleaner
    transitions = plasma.atomic_data.macro_atom_data[
        plasma.atomic_data.macro_atom_data.transition_type == -1
    ].copy()
    transitions_index = transitions.set_index(
        ["atomic_number", "ion_number", "source_level_number"]
    ).index.copy()
    tmp = pd.DataFrame(
        transition_probabilities[
            (atomic_data.macro_atom_data.transition_type == -1).values
        ]
    )
    q_ul = tmp.set_index(transitions_index)
    t = simulation_state.time_explosion.value
    lines = atomic_data.lines.set_index("line_id")
    wave = lines.wavelength_cm.loc[
        transitions.transition_line_id
    ].values.reshape(-1, 1)
    if transport.line_interaction_type == "macroatom":
        e_dot_u = C_frame.loc[e_dot_u.index]
    att_S_ul = wave * (q_ul * e_dot_u) * t / (4 * np.pi)

    result = pd.DataFrame(
        att_S_ul.values,
        index=transitions.transition_line_id.values,
        columns=columns,
    )
    att_S_ul = result.loc[lines.index.values].values

    # Jredlu should already by in the correct order, i.e. by wavelength of
    # the transition l->u (similar to Jbluelu)
    Jredlu = Jbluelu * np.exp(-tau_sobolevs) + att_S_ul
    if interpolate_shells > 0:
        (
            att_S_ul,
            Jredlu,
            Jbluelu,
            e_dot_u,
        ) = interpolate_integrator_quantities(
            att_S_ul, Jredlu, Jbluelu, e_dot_u,
            interpolate_shells,
            transport, simulation_state, opacity_state, plasma
        )
    else:
        transport.r_inner_i = (
            montecarlo_transport_state.geometry_state.r_inner
        )
        transport.r_outer_i = (
            montecarlo_transport_state.geometry_state.r_outer
        )
        transport.tau_sobolevs_integ = opacity_state.tau_sobolev
        transport.electron_densities_integ = (
            opacity_state.electron_density
        )

    return att_S_ul, Jredlu, Jbluelu, e_dot_u

def interpolate_integrator_quantities(
    att_S_ul, Jredlu, Jbluelu, e_dot_u,
    interpolate_shells,
    transport, simulation_state, opacity_state, plasma
):

    mct_state = transport.transport_state
    
    nshells = interpolate_shells
    r_middle = (
        mct_state.geometry_state.r_inner + mct_state.geometry_state.r_outer
    ) / 2.0

    r_integ = np.linspace(
        mct_state.geometry_state.r_inner[0],
        mct_state.geometry_state.r_outer[-1],
        nshells,
    )
    transport.r_inner_i = r_integ[:-1]
    transport.r_outer_i = r_integ[1:]

    r_middle_integ = (r_integ[:-1] + r_integ[1:]) / 2.0

    transport.electron_densities_integ = interp1d(
        r_middle,
        plasma.electron_densities.iloc[
            simulation_state.geometry.v_inner_boundary_index : simulation_state.geometry.v_outer_boundary_index
        ],
        fill_value="extrapolate",
        kind="nearest",
    )(r_middle_integ)
    # Assume tau_sobolevs to be constant within a shell
    # (as in the MC simulation)
    transport.tau_sobolevs_integ = interp1d(
        r_middle,
        opacity_state.tau_sobolev[
            :,
            simulation_state.geometry.v_inner_boundary_index : simulation_state.geometry.v_outer_boundary_index,
        ],
        fill_value="extrapolate",
        kind="nearest",
    )(r_middle_integ)
    att_S_ul = interp1d(r_middle, att_S_ul, fill_value="extrapolate")(
        r_middle_integ
    )
    Jredlu = interp1d(r_middle, Jredlu, fill_value="extrapolate")(
        r_middle_integ
    )
    Jbluelu = interp1d(r_middle, Jbluelu, fill_value="extrapolate")(
        r_middle_integ
    )
    e_dot_u = interp1d(r_middle, e_dot_u, fill_value="extrapolate")(
        r_middle_integ
    )

    # Set negative values from the extrapolation to zero
    att_S_ul = att_S_ul.clip(0.0)
    Jbluelu = Jbluelu.clip(0.0)
    Jredlu = Jredlu.clip(0.0)
    e_dot_u = e_dot_u.clip(0.0)
    return att_S_ul, Jredlu, Jbluelu, e_dot_u
