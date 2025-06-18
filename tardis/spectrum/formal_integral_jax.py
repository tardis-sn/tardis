import jax
from jax import jit
import jax.numpy as jnp
jax.config.update("jax_enable_x64", True)
from functools import partial

import timeit

import numpy as np
import math 

C_INV = 3.33564e-11
PI = np.pi
KB_CGS = 1.3806488e-16
H_CGS = 6.62606957e-27

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

@partial(jit, static_argnames='inv_t')
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
    return jnp.where(r > p,
              jnp.sqrt(r * r - p * p) * C_INV * inv_t,
              0.0)

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
    inv_t = 1./time_explosion

    def in_photosphere(arrs):
        z0, shell_id = arrs

        idx = jnp.arange(0, N)
        z0 = z0.at[idx].set(1-calculate_z_jax(r, p, inv_t))

        shell_id = shell_id.at[idx].set(idx)
        return z0, shell_id
    
    def out_photosphere_old(arrs):

        # TODO: not getting the offset correctly
        z0, shell_id = arrs

        z = calculate_z_jax(r, p, inv_t)
        
        # otherwise get the edges of the shell
        idx_low = jnp.arange(N - 1, -1, -1) 
        z0 = z0.at[idx_low].set(1+z)
        shell_id = shell_id.at[idx_low].set(jnp.arange(0, N))

        # offset is the first i that z is non-zero
        offset = N

        idx_up = jnp.arange(N, N+N)
        z0 = z0.at[idx_up].set(1-z)
        shell_id = shell_id.at[idx_up].set(jnp.arange(0, N))

        # replace zeroes if z was zero
        idx_zeroes = jnp.where(z == 0, size=N)[0] 
        z0 = z0.at[idx_low[idx_zeroes]].set(0)
        z0 = z0.at[idx_up[idx_zeroes]].set(0)
        shell_id = shell_id.at[idx_low[idx_zeroes]].set(0)
        shell_id = shell_id.at[idx_up[idx_zeroes]].set(0)
        
        return z0, shell_id

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
        z0, shell_id, _ = jax.lax.fori_loop(0, N, loop_out, (z0, shell_id, N))
        return z0, shell_id
    
    cond = p <= r_inner[0]
    # print(cond)
    return jax.lax.cond(cond, in_photosphere, out_photosphere, (z0, shell_id))

# compute for all p
calc_z_jax_jitvmap = jax.jit(jax.vmap(populate_z_jax, in_axes=[0, None, None, None, None]), static_argnames=['time_explosion', 'size_shell'])

def intensity_black_body_jax(nu, temperature):
    beta_rad = 1 / (KB_CGS * temperature)
    coefficient = 2 * H_CGS * C_INV * C_INV

    return jnp.where(nu == 0, jnp.nan, coefficient * nu * nu * nu / (jnp.exp(H_CGS * nu * beta_rad) - 1))

# in_axes done like this indicates that the 2nd arg is fixed in the vmaping
ibb_vmap = jax.jit(jax.vmap(intensity_black_body_jax, in_axes=[0, None]))

@partial(jit, static_argnames=['imin', 'imax'])
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

    cond = (x_insert > x[imin]) or (x_insert < x[imax])
    return jnp.where(cond, -1, len(x) - 1 - jnp.searchsorted(x[::-1], x_insert, side="right"))

@partial(jit, static_argnames='nu')
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

    def binary_search(_):
        result = reverse_binary_search_jax(nu, nu_insert, imin, imax)
        return result + 1

    return jax.lax.cond(nu_insert > nu[imin],
                0.,
                jax.lax.cond(nu_insert < nu[imax],
                                imax + 1,
                                binary_search))