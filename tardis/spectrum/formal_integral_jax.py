import jax
from jax import jit
import jax.numpy as jnp
import jax.debug as jdb

jax.config.update("jax_enable_x64", True)
from functools import partial

import timeit

import numpy as np
import math
from scipy.interpolate import interp1d

from tardis.transport.montecarlo.configuration.constants import SIGMA_THOMSON

C_INV = 3.33564e-11
PI = np.pi
KB_CGS = 1.3806488e-16
H_CGS = 6.62606957e-27

@partial(jit, static_argnames="N")
def calculate_p_values_jax(Rmax, N):
    """
    Calculates N p values

    Parameters
    ----------
    Rmax : float64
        maximum radius
    N : int64
        number of impact parameters

    Returns
    -------
    jnp.array float64
    """
    return jnp.arange(N).astype(jnp.float64) * Rmax / (N - 1)


@partial(jit, static_argnames="inv_t")
def calculate_z_jax(r, p, inv_t):
    """
    Calculate distances to p lines

    Calculate half of the length of the p-line inside a shell
    of radius r in terms of unit length (c * t_exp).
    If shell and p-line do not intersect, return 0.

    Parameters
    ----------
    r : jnp.array float64
        outer radii of the shells
    p : jnp.array float64
        distances of the p-lines to the center of the supernova
    inv_t : float64
        inverse time_explosion

    Returns
    -------
    float64
    """
    return jnp.where(r > p, jnp.sqrt(r * r - p * p) * C_INV * inv_t, 0.0)


def populate_z_jax(p, r_outer, r_inner, time_explosion, size_shell):
    """
    Calculates the intersection points of the p-line with each shell

    Parameters
    ----------
    p : jnp.array
        distances of the p-lines to the center of the supernova
    r_outer : jnp.array
        outer radii of the shells
    r_inner :jnp.array
        inner radii of the shells
    time_explosion : float64
        time of the explosion in seconds
    size_shell : _type_
        length of the shell array

    Returns
    -------
    jnp.array float64
    """

    z0 = jnp.zeros(2 * size_shell, dtype=jnp.float64)
    shell_id = np.zeros(2 * size_shell, dtype=np.int64)

    r = r_outer
    N = len(r_inner)
    inv_t = 1.0 / time_explosion

    def in_photosphere(arrs):
        z0, shell_id = arrs

        idx = jnp.arange(0, N)
        z0 = z0.at[idx].set(1 - calculate_z_jax(r, p, inv_t))

        shell_id = shell_id.at[idx].set(idx)
        return z0, shell_id, N

    # compute interactions if outside of the photosphere

    # loop over the radii, update z0 if z is not zero
    def loop_out(i, state):
        ri = r[i]
        z = calculate_z_jax(r[i], p, inv_t)

        def do_update(state):
            z0, shell_id, offset = state
            offset = jnp.where(offset == N, i, offset)
            i_low = N - i - 1
            i_up = N + i - 2 * offset

            z0 = z0.at[i_low].set(1 + z)
            z0 = z0.at[i_up].set(1 - z)
            shell_id = shell_id.at[i_low].set(i)
            shell_id = shell_id.at[i_up].set(i)

            return z0, shell_id, offset

        def do_nothing(state):
            return state

        state = jax.lax.cond(z == 0, do_nothing, do_update, state)
        return state

    def out_photosphere(arrs):
        z0, shell_id = arrs
        z0, shell_id, offset = jax.lax.fori_loop(
            0, N, loop_out, (z0, shell_id, N)
        )
        return z0, shell_id, 2 * (N - offset)

    cond = p <= r_inner[0]
    return jax.lax.cond(cond, in_photosphere, out_photosphere, (z0, shell_id))


# jit compile and vectorize over impact parameters
calc_z_jax_jitvmap = jax.jit(
    jax.vmap(populate_z_jax, in_axes=[0, None, None, None, None]),
    static_argnames=["time_explosion", "size_shell"],
)


@jit
def intensity_black_body_jax(nu, temperature):
    beta_rad = 1 / (KB_CGS * temperature)
    coefficient = 2 * H_CGS * C_INV * C_INV

    return jnp.where(
        nu == 0,
        jnp.nan,
        coefficient * nu * nu * nu / (jnp.exp(H_CGS * nu * beta_rad) - 1),
    )


@partial(jit, static_argnames=["x", "imin", "imax"])
def reverse_binary_search_jax(x, x_insert, imin, imax):
    """
    Find indicies where elements should be inserted
    to maintain order in an inversely sorted float
    array.

    Find the indices into a sorted array a such that,
    if the corresponding elements in v were inserted
    before the indices on the right, the order of a
    would be preserved.

    Parameters
    ----------
    nu : jnp.array(jnp.float64, 1d)
        line frequencies
    nu_insert : int
        value of nu key
    imin: int
        Lower bound
    imax: int
        Upper bound

    Returns
    -------
    int
        Location of insertion
    """

    # cond = (x_insert > x[imin]) & (x_insert < x[imax])
    # res = jnp.where(
    #     cond, -1, len(x) - 1 - jnp.searchsorted(x[::-1], x_insert, side="right")
    # )
    return len(x) - 1 - jnp.searchsorted(x[::-1], x_insert, side="right")
    

    # TODO: warn about lack of checking
    # return len(x) - 1 - jnp.searchsorted(x[::-1], x_insert, side="right")

# @partial(jax.vmap, in_axes=(None, 0))
@jit
def line_search_jax(nu, nu_insert):
    """
    Insert a value in to an array of line frequencies

    Parameters
    ----------
    nu : jnp.array(jnp.float64)
        line frequencies
    nu_insert : int
        value of nu key

    Returns
    -------
    int
        index of the next line to the red.
        If the key value is redder than
        the reddest line returns number_of_lines.
    """

    imin = 0
    imax = len(nu) - 1

    # TODO: check if this gets compiled
    # TODO: check the types of these things
    def binary_search():
        result = reverse_binary_search_jax(nu, nu_insert, imin, imax)
        return jnp.int64(result + 1)

    def no_search():
        return jnp.int64(imax + 1)

    return jax.lax.cond(
        nu_insert > nu[imin],
        lambda: jnp.int64(0), # must be callable
        lambda: jax.lax.cond(nu_insert < nu[imax], no_search, binary_search),
    )


def init_Inup(nu, p, z, iT, Rph):
    val = jnp.where(
        (p <= Rph) & (p != 0.0), intensity_black_body_jax(nu * z, iT), 0.0
    )
    return val
init_Inup_jax = jax.jit(
    jax.vmap(
        jax.vmap(init_Inup, in_axes=(None, 0, 0, None, None)),
        in_axes=(0, None, None, None, None)
    )
)

# inner for loop:
def calc_Inup(
    Inup,
    nu,
    p,
    z,
    size_z,
    shell_id,
    line_list_nu,
    time_explosion,
    electron_density,
    exp_tau,
    att_S_ul,
    Jred_lu,
    Jblue_lu,
):
    
    size_line = len(line_list_nu)
    nu_start = nu * z[0]
    pline = line_search_jax(
        line_list_nu, nu_start
    )  # ensure this computes all plines
    offset = shell_id[0] * size_line
    pline_offset = pline + offset

    nu_ends = nu * z[1:]
    nu_ends_idxs = (
        size_line
        - 1
        - jnp.searchsorted(line_list_nu[::-1], nu_ends, side="right")
    )

    # inner loop
    def calc_escat(i, escat_params):
        # first contribution to integration
        def first_contrib(cond_state):
            escat_contrib, zdiff, escat_op, Inup, pline_offset = cond_state
            escat_contrib += zdiff * escat_op * (Jblue_lu[pline_offset] - Inup)
            return escat_contrib, 0

        # Account for e-scattering, c.f. Eqs 27, 28 in Lucy 1999
        def otherwise(cond_state):
            escat_contrib, zdiff, escat_op, Inup, pline_offset = cond_state

            Jkkp = 0.5 * (Jred_lu[pline_offset - 1] + Jblue_lu[pline_offset])
            escat_contrib += (zdiff) * escat_op * (Jkkp - Inup)
            return escat_contrib, 0

        Inup, pline, pline_offset, escat_contrib, escat_op, zstart, first = escat_params
        zend = time_explosion / C_INV * (1.0 - line_list_nu[pline] / nu)
        zdiff = zend - zstart
        escat_contrib, first = jax.lax.cond(
            first == 1,
            first_contrib,
            otherwise,
            (0, zdiff, escat_op, Inup, pline_offset),
        )
    
        Inup += escat_contrib
        Inup *= exp_tau[pline_offset]
        Inup += att_S_ul[pline_offset]

        # increment indexed and reset zstart
        return (Inup, pline + 1, pline_offset + 1, escat_contrib, escat_op, zend, first)

    def loop_interactions(i, loop_params):
        Inup, pline, pline_offset, escat_contrib, zstart = loop_params
        escat_op = electron_density[shell_id[i]] * SIGMA_THOMSON
        nu_end = nu_ends[i]
        nu_end_idx = nu_ends_idxs[i]
        
        # escat loop
        Inup, pline, pline_offset, _, _, zstart, _ = jax.lax.fori_loop(
            0,
            jnp.maximum(nu_end_idx - pline, 0),
            calc_escat,
            (Inup, pline, pline_offset, escat_contrib, escat_op, zstart, 1),
        )

        # calculate e-scattering optical depth to grid cell boundary
        Jkkp = 0.5 * (Jred_lu[pline_offset - 1] + Jblue_lu[pline_offset])
        zend = time_explosion / C_INV * (1.0 - nu_end / nu)
        escat_contrib += (zend - zstart) * escat_op * (Jkkp - Inup)
        dir = (shell_id[i + 1] - shell_id[i]) * size_line # TODO: replace since flattened need to move by sizeline, but if not flat then just diff axis

        return (Inup, pline + dir, pline_offset + dir, escat_contrib, zend)

    # run the inner loop for i in range(size_z - 1):
    zstart = time_explosion / C_INV * (1.0 - z[0])
    Inup, _, _, _, _ = jax.lax.fori_loop(
        0,
        size_z - 1,
        loop_interactions,
        (Inup, pline, pline_offset, 0.0, zstart),
    )

    Inup *= p
    return Inup

# TODO: determine how static variables work
    # it complained about the jax arrays being unhashable but I want it to know that the first time it sees those arrays they will always be that
calc_Inup_vmapped = jax.vmap(
    jax.vmap(
        calc_Inup,
        in_axes=(0, None, 0, 0, 0, 0, None, None, None, None, None, None, None)
        #        ^Inup[nu,p] ^nu ^p ^z[p,:] ^size_z ^shell_id[p,:] ^line_list_nu ^time_explosion ^electron_density ^exp_tau ^att_S_ul ^Jred_lu ^Jblue_lu
    ),
    in_axes=(0, 0, None, None, None, None, None, None, None, None, None, None,None)
    #        ^Inup[nu,:] ^nu ^p ^z ^size_z ^shell_id ^line_list_nu ^time_explosion ^electron_density ^exp_tau ^att_S_ul ^Jred_lu ^Jblue_lu
)
calc_Inup_jax = jax.jit(calc_Inup_vmapped)


@partial(jit, static_argnames=["interpolate_shells", "v_inner_idx", "v_outer_idx"])
def interpolate_integrator_quantities_jax(mct_r_inner, mct_r_outer, v_inner_idx, v_outer_idx, 
                                      tau_sobolev, electron_densities, 
                                      interpolate_shells, att_S_ul, Jredlu, Jbluelu, e_dot_u):
    
    r_middle = (
        mct_r_inner + mct_r_outer
    ) / 2.0 # is this a different geometry?

    r_integ = jnp.linspace(
        mct_r_inner[0],
        mct_r_outer[-1],
        interpolate_shells,
    )

    r_inner_i = r_integ[:-1]
    r_outer_i = r_integ[1:]

    r_middle_integ = (r_integ[:-1] + r_integ[1:]) / 2.0

    # interp values
    electron_densities_integ = jax.scipy.interpolate.RegularGridInterpolator(
            (r_middle, ),
            electron_densities[v_inner_idx:v_outer_idx],
            fill_value=None, # TODO: determine fill
            method="nearest", 
        )(r_middle_integ)
    
    # since these array are of shape (nu, shells), they have to be transposed such that the shells are first
    tau_sobolevs_integ = jax.scipy.interpolate.RegularGridInterpolator(
            (r_middle,),
            tau_sobolev[:, v_inner_idx:v_outer_idx].T,
            method="nearest",
            fill_value=None  # for extrapolation    
        )(r_middle_integ.reshape(-1, 1)).T
    att_S_ul = jax.scipy.interpolate.RegularGridInterpolator(
            (r_middle,), att_S_ul.T, fill_value=None
        )(r_middle_integ.reshape(-1, 1)).T
    Jredlu = jax.scipy.interpolate.RegularGridInterpolator(
            (r_middle,), Jredlu.T, fill_value=None
        )(r_middle_integ.reshape(-1, 1)).T
    Jbluelu = jax.scipy.interpolate.RegularGridInterpolator(
            (r_middle,), Jbluelu.T, fill_value=None
        )(r_middle_integ.reshape(-1, 1)).T
    e_dot_u = jax.scipy.interpolate.RegularGridInterpolator(
            (r_middle,), e_dot_u.T, fill_value=None
        )(r_middle_integ.reshape(-1, 1)).T

    # Set negative values from the extrapolation to zero
    att_S_ul = jnp.clip(att_S_ul, min=0.0)
    Jbluelu = jnp.clip(Jbluelu, min=0.0)
    Jredlu = jnp.clip(Jredlu, min=0.0)
    e_dot_u = jnp.clip(e_dot_u, min=0.0)

    return r_inner_i, r_outer_i, electron_densities_integ, tau_sobolevs_integ, att_S_ul, Jredlu, Jbluelu, e_dot_u


# TODO: probably cant pass in the objects if I want to jjax
# TODO: state assumption about only macroatom and not downbranch - reduces conditionals
def make_source_function(simulation_state, transport, opacity_state, atomic_data, levels_index):

    local_slice = slice(
            simulation_state.geometry.v_inner_boundary_index,
            simulation_state.geometry.v_outer_boundary_index,
    )
    montecarlo_transport_state = transport.transport_state
    transition_probabilities = opacity_state.transition_probabilities[
            :, local_slice
    ]
    tau_sobolevs = opacity_state.tau_sobolev[:, local_slice]

    columns = range(simulation_state.no_of_shells)

    macro_ref = atomic_data.macro_atom_reference
    macro_data = atomic_data.macro_atom_data

    no_lvls = len(levels_index)
    no_shells = len(simulation_state.dilution_factor)


    return 



def formal_integral(
    r_inner,
    r_outer,
    time_explosion,
    tau_sobolev,
    line_list_nu,
    iT,
    inu,
    att_S_ul,
    Jred_lu,
    Jblue_lu,
    electron_density,
    N,
):
    """
    Returns
    -------
    L : float64 array
        integrated luminosities
    I_nu_p : float64 2D array
        intensities at each p-ray multiplied by p
        frequency x p-ray grid
    """
    # get params
    size_line, size_shell = tau_sobolev.shape
    exp_tau = jnp.exp(-tau_sobolev.T.ravel()) # TODO: figure this out (10000 nus but 29000 in first axis)

    # compute impact parameters, p
    ps = calculate_p_values_jax(r_outer[-1], N)

    # compute interaction points, z
    zs, shellids, size_zs = calc_z_jax_jitvmap(
        ps, r_outer, r_inner, time_explosion, size_shell
    )

    # init Inup with the photopshere's black body
    Inup = init_Inup_jax(inu, ps, zs[:, 0], iT, r_inner[0])

    Inup_i = calc_Inup_jax(
        Inup[:, 1:],
        inu,
        ps[1:],
        zs[1:],
        size_zs[1:],
        shellids[1:], # #TODO: determine if correct (unsure)
        line_list_nu,
        time_explosion,
        electron_density,
        exp_tau,
        att_S_ul.ravel(order="F"),
        Jred_lu.ravel(order="F"),
        Jblue_lu.ravel(order="F"),
    )

    return Inup_i
