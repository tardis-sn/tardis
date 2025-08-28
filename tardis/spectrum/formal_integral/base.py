import numpy as np
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


def interpolate_integrator_quantities(
    source_function_state,
    interpolate_shells,
    simulation_state,
    transport,
    opacity_state,
    electron_densities,
):
    """Interpolate the integrator quantities to interpolate_shells.

    Parameters
    ----------
    source_function_state: SourceFunctionState
        Data class to hold the computed source function values
    interpolate_shells : int
        number of shells to interpolate to
    simulation_state : tardis.model.SimulationState
    transport : tardis.transport.montecarlo.MonteCarloTransportSolver
    opacity_state : OpacityStateNumba
    electron_densities : pd.DataFrame

    Returns
    -------
    tuple
        Interpolated values of att_S_ul, Jred_lu, Jbluelu, and e_dot_u
    """

    mct_state = transport.transport_state

    att_S_ul, Jred_lu, Jblue_lu, e_dot_u = (
            source_function_state.att_S_ul,
            source_function_state.Jred_lu,
            source_function_state.Jblue_lu,
            source_function_state.e_dot_u,
        )

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
        electron_densities.iloc[
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
    Jred_lu = interp1d(r_middle, Jred_lu, fill_value="extrapolate")(
        r_middle_integ
    )
    Jblue_lu = interp1d(r_middle, Jblue_lu, fill_value="extrapolate")(
        r_middle_integ
    )
    e_dot_u = interp1d(r_middle, e_dot_u, fill_value="extrapolate")(
        r_middle_integ
    )

    # Set negative values from the extrapolation to zero
    att_S_ul = att_S_ul.clip(0.0)
    Jblue_lu = Jblue_lu.clip(0.0)
    Jred_lu = Jred_lu.clip(0.0)
    e_dot_u = e_dot_u.clip(0.0)
    return att_S_ul, Jred_lu, Jblue_lu, e_dot_u


import numpy as np
from numba import njit, prange

from tardis.spectrum.formal_integral.base import (
    C_INV,
    BoundsError,
    intensity_black_body,
    calculate_p_values,
)
from tardis.transport.montecarlo import njit_dict, njit_dict_no_parallel
from tardis.transport.montecarlo.configuration.constants import SIGMA_THOMSON

@njit(**njit_dict_no_parallel)
def calculate_z(r, p, inv_t):
    """Calculate distance to p line

    Calculate half of the length of the p-line inside a shell
    of radius r in terms of unit length (c * t_exp).
    If shell and p-line do not intersect, return 0.

    Inputs:
        :r: (double) radius of the shell
        :p: (double) distance of the p-line to the center of the supernova
        :inv_t: (double) inverse time_explosion is needed to norm to unit-length
    """
    if r > p:
        return np.sqrt(r * r - p * p) * C_INV * inv_t
    else:
        return 0


@njit(**njit_dict_no_parallel)
def populate_z(geometry, time_explosion, p, oz, oshell_id):
    """Calculate p line intersections

    This function calculates the intersection points of the p-line with
    each shell

    Inputs:
        :p: (double) distance of the integration line to the center
        :oz: (array of doubles) will be set with z values. the array is truncated
                    by the value `1`.
        :oshell_id: (int64) will be set with the corresponding shell_ids
    """
    # abbreviations
    r = geometry.r_outer
    N = len(geometry.r_inner)  # check
    inv_t = 1 / time_explosion
    z = 0
    offset = N

    if p <= geometry.r_inner[0]:
        # intersect the photosphere
        for i in range(N):
            oz[i] = 1 - calculate_z(r[i], p, inv_t)
            oshell_id[i] = i
        return N
    else:
        # no intersection with photosphere
        # that means we intersect each shell twice
        for i in range(N):
            z = calculate_z(r[i], p, inv_t)
            if z == 0:
                continue
            if offset == N:
                offset = i
            # calculate the index in the resulting array
            i_low = N - i - 1  # the far intersection with the shell
            i_up = N + i - 2 * offset  # the nearer intersection with the shell

            # setting the arrays; check return them?
            oz[i_low] = 1 + z
            oshell_id[i_low] = i
            oz[i_up] = 1 - z
            oshell_id[i_up] = i
        return 2 * (N - offset)


@njit(**njit_dict_no_parallel)
def reverse_binary_search(x, x_insert, imin, imax):
    """Look for a place to insert a value in an inversely sorted float array.

    Inputs:
        :x: (array) an inversely (largest to lowest) sorted float array
        :x_insert: (value) a value to insert
        :imin: (int) lower bound
        :imax: (int) upper bound

    Outputs:
        index of the next boundary to the left
    """
    # ret_val = TARDIS_ERROR_OK # check
    if (x_insert > x[imin]) or (x_insert < x[imax]):
        raise BoundsError  # check
    return len(x) - 1 - np.searchsorted(x[::-1], x_insert, side="right")


@njit(**njit_dict_no_parallel)
def line_search(nu, nu_insert, number_of_lines):
    """
    Insert a value in to an array of line frequencies

    Inputs:
        :nu: (array) line frequencies
        :nu_insert: (int) value of nu key
        :number_of_lines: (int) number of lines in the line list

    Outputs:
        index of the next line ot the red.
                If the key value is redder
                 than the reddest line returns number_of_lines.
    """
    # TODO: fix the TARDIS_ERROR_OK
    # tardis_error_t ret_val = TARDIS_ERROR_OK # check
    imin = 0
    imax = number_of_lines - 1
    if nu_insert > nu[imin]:
        result = imin
    elif nu_insert < nu[imax]:
        result = imax + 1
    else:
        result = reverse_binary_search(nu, nu_insert, imin, imax)
        result = result + 1
    return result


@njit(**njit_dict_no_parallel)
def trapezoid_integration(array, h):
    """in the future, let's just replace
    this with the numpy trapz
    since it is numba compatable
    """
    return np.trapz(array, dx=h)


# numba jit
calculate_p_values = njit(calculate_p_values, **njit_dict_no_parallel)
intensity_black_body = njit(intensity_black_body, **njit_dict_no_parallel)


# initialize I_nu_p with the blackbody intensities (as necessary)
@njit(**njit_dict)
def init_Inup_numba(inu, ps, iT, Rph, N, size_shell, geometry, time_explosion):
    inu_size = len(inu)

    I_nu_p = np.zeros((inu_size, N), dtype=np.float64)
    zs = np.zeros((N, 2 * size_shell), dtype=np.float64)
    shell_ids = np.zeros((N, 2 * size_shell), dtype=np.int64)
    size_zs = np.zeros(N, dtype=np.int64)

    # loop over nu and p
    for nu_idx in prange(inu_size):
        I_nu = I_nu_p[nu_idx]
        nu = inu[nu_idx]
        for p_idx in range(1, N):
            p = ps[p_idx]

            # get zs
            size_z = populate_z(
                geometry, time_explosion, p, zs[p_idx], shell_ids[p_idx]
            )
            size_zs[p_idx] = size_z

            # if inside the photosphere, set to black body intensity
            # otherwise zero
            if p <= Rph:
                I_nu[p_idx] = intensity_black_body(nu * zs[p_idx][0], iT)
            else:
                I_nu[p_idx] = 0

    return I_nu_p, zs, shell_ids, size_zs

@njit(**njit_dict)
def setup_formal_integral_inputs(
        inu,
        iT,
        N, 
        geometry, 
        time_explosion,
        line_list_nu,
        tau_sobolev, 
        electron_densities
    ):
    """
    Prepare all of the arrays and values needed for the loops inside the formal integral

    Needed parameters at the end:
        nu, ps, Inup
        att_S_ul, Jred_lu, Jblue_lu, exp_tau
        escat_op
        line_list
        shell_ids, size_lines,
        direction (?) can also be precomputed
    """

    inu_size = len(inu)
    size_line, size_shell = tau_sobolev.shape
    exp_tau = np.exp(-tau_sobolev.T.ravel())
    R_max = geometry.r_outer[size_shell - 1]
    Rph = geometry.r_inner[0]

    ps = np.zeros(N, dtype=np.float64)
    ps[::] = calculate_p_values(R_max, N)

    I_nu_p = np.zeros((inu_size, N), dtype=np.float64)
    zs = np.zeros((N, 2 * size_shell), dtype=np.float64)
    shell_ids = np.zeros((N, 2 * size_shell), dtype=np.int64)
    size_zs = np.zeros(N, dtype=np.int64)

    # nu start and nu_start idx not needed other than to compute internal values
    lines_idx = np.zeros((inu_size, N), dtype=np.int64)
    lines_idx_offset = np.zeros((inu_size, N), dtype=np.int64)
    
    zstarts = np.zeros(N, dtype=np.float64)

    nu_ends = np.zeros((inu_size, N, 2 * size_shell-1), dtype=np.float64)
    nu_ends_idxs = np.zeros((inu_size, N, 2 * size_shell-1), dtype=np.int64)

    escat_ops = np.zeros((inu_size, N, 2 * size_shell), dtype=np.float64)
    directions = np.zeros((inu_size, N, 2 * size_shell), dtype=np.int64)

    # loop over nu and p
    for nu_idx in prange(inu_size):
        I_nu = I_nu_p[nu_idx]
        nu = inu[nu_idx]
        for p_idx in range(1, N):
            p = ps[p_idx]

            # get zs
            size_z = populate_z(
                geometry, time_explosion, p, zs[p_idx], shell_ids[p_idx]
            )
            size_zs[p_idx] = size_z
            z = zs[p_idx]
            shell_id = shell_ids[p_idx]

            # if inside the photosphere, set to black body intensity
            # otherwise zero
            if p <= Rph:
                I_nu[p_idx] = intensity_black_body(nu * z[0], iT)
            else:
                I_nu[p_idx] = 0

            nu_start = nu * z[0]

            idx_nu_start = line_search(line_list_nu, nu_start, size_line)
            lines_idx[nu_idx, p_idx] = int(idx_nu_start)
            lines_idx_offset[nu_idx, p_idx] = int(idx_nu_start + (shell_id[0] * size_line))
            
            zstarts[p_idx] = time_explosion / C_INV * (1.0 - z[0]) 

            nu_ends[nu_idx, p_idx] = nu * z[1:]
            nu_ends_idxs[nu_idx, p_idx] = size_line - np.searchsorted(
                line_list_nu[::-1], nu_ends[nu_idx, p_idx], side="right"
            )
            
            for i in range(size_z - 1):
                escat_ops[nu_idx, p_idx, i] = (
                    electron_densities[int(shell_ids[p_idx][i])] * SIGMA_THOMSON
                )
                # compute nu_end_idx - line_idx
                directions[nu_idx, p_idx, i] = int((shell_ids[p_idx][i + 1] - shell_ids[p_idx][i]) * size_line)


    return I_nu_p, ps, size_zs, lines_idx, lines_idx_offset, nu_ends, nu_ends_idxs, zstarts, escat_ops, directions, exp_tau