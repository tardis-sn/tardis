from functools import partial
import timeit

import jax
from jax import jit
import jax.numpy as jnp
import jax.debug as jdb

jax.config.update("jax_enable_x64", True)

import numpy as np
from astropy import units as u

from tardis.transport.montecarlo.configuration.constants import SIGMA_THOMSON
from tardis.spectrum.formal_integral.base import C_INV, KB_CGS, H_CGS

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


def populate_z(p, r_inner, r_outer, time_explosion, size_shell):
    """
    Calculates the intersection points of the p-line with each shell

    Parameters
    ----------
    p : jnp.array float64
        distances of the p-lines to the center of the supernova
    r_outer : jnp.array float64
        outer radii of the shells
    r_inner :jnp.array float64
        inner radii of the shells
    time_explosion : float64
        time of the explosion in seconds
    size_shell : int64
        length of the shell array

    Returns
    -------
    z: jnp.array float64
        the intersection points
    shell_id: jnp.array int64
        ids of the shell intersected
    size_z: int64
        number of interaction points
    """

    z0 = jnp.zeros(2 * size_shell, dtype=jnp.float64)
    shell_id = np.zeros(2 * size_shell, dtype=np.int64)

    r = r_outer
    N = len(r_inner)
    inv_t = 1.0 / time_explosion

    # compute interactions if inside the photosphere
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
        z = calculate_z_jax(ri, p, inv_t)

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


# jit compile and vectorize over impact parameters p
populate_z_jax = jax.jit(
    jax.vmap(populate_z, in_axes=[0, None, None, None, None]),
    static_argnames=["time_explosion", "size_shell"],
)


@jit
def intensity_black_body_jax(nu, temperature):
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

    beta_rad = 1 / (KB_CGS * temperature)
    coefficient = 2 * H_CGS * C_INV * C_INV

    return jnp.where(
        nu == 0,
        jnp.nan,
        coefficient * nu * nu * nu / (jnp.exp(H_CGS * nu * beta_rad) - 1),
    )


@partial(jit, static_argnames=["imin", "imax"])
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
    # TODO: warn about lack of checking
    return len(x) - 1 - jnp.searchsorted(x[::-1], x_insert, side="right")


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
        lambda: jnp.int64(0),  # must be callable
        lambda: jax.lax.cond(nu_insert < nu[imax], no_search, binary_search),
    )


def init_Inup(nu, p, z, iT, Rph):
    """
    Creates the initial intensity per nu and p array

    Parameters
    ----------
    nu : float64
        frequency
    p : float64
        impact parameter
    z : array float64
        interaction points
    iT : float64
        interpolated temperature
    Rph : float64
        radius of the photosphere

    Returns
    -------
    float64
        initial intensity per nu and p
    """

    val = jnp.where(
        (p <= Rph) & (p != 0.0), intensity_black_body_jax(nu * z, iT), 0.0
    )
    return val

# jit compile and vectorize over nu and p
init_Inup_jax = jax.jit(
    jax.vmap(
        jax.vmap(init_Inup, in_axes=(None, 0, 0, None, None)),
        in_axes=(0, None, None, None, None),
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
    """
    Updates the initial intensity per nu, p, and interaction

    Parameters
    ----------
    Inup : float64
        initial intensity per nu and p
    nu : float64
        frequency
    p : float64
        impact parameter
    z : float64
        intersection points
    size_z : int
        number of interaction points
    shell_id : int
        ids of the shell intersected
    line_list_nu : float64 array
        line frequencies
    time_explosion : float64
        time of the explosion in seconds
    electron_density : float64 array
        electron densities in each shell
    exp_tau : float64 array
        exponential of the Sobolev optical depth
    att_S_ul : float64 array
        attenuated source function for each line in each shell
    Jred_lu : float64 array
        J estimator from red end of the line from lower to upper level
    Jblue_lu : float64 array
        J estimator from blue end of the line from lower to upper level

    Returns
    -------
    Inup : float64
    """

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

        Inup, pline, pline_offset, escat_contrib, escat_op, zstart, first = (
            escat_params
        )
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
        return (
            Inup,
            pline + 1,
            pline_offset + 1,
            escat_contrib,
            escat_op,
            zend,
            first,
        )

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
        dir = (
            (shell_id[i + 1] - shell_id[i]) * size_line
        )  # TODO: replace since flattened need to move by sizeline, but if not flat then just diff axis

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
calc_Inup_jax = jax.jit(
    jax.vmap(  # vectorize over nu
        jax.vmap(  # vectorize over p
            calc_Inup,
            in_axes=(
                0,
                None,
                0,
                0,
                0,
                0,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            ),
        ),
        in_axes=(
            0,
            0,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        ),
    )
)
# p_jax = jax.jit(calc_Inup_vmapped)

def formal_integral_jax(
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
    Computes the formal integral

    Parameters
    ----------
    r_inner : float64 array
        inner radii of the shells
    r_outer : float64 array
        outer radii of the shells
    time_explosion : float64
        time of the explosion in seconds
    tau_sobolev : float64 2D array
        Sobolev Optical depth for each line in each shell
    line_list_nu : float64 array
        line frequencies
    iT : float64
        interpolated temperature
    inu : float64 array
        interpolated frequencies
    att_S_ul : float64 2D array
        attentuated source function
    Jred_lu : float64 2D array
        J estimator from red end of the line from lower to upper level
    Jblue_lu : float64 2D array
        J estimator from blue end of the line from lower to upper level
    electron_density : float64 1D array
        electron densities in each shell
    N : int
        number of impact parameters


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
    exp_tau = jnp.exp(
        -tau_sobolev.T.ravel()
    )  # TODO: figure this out (10000 nus but 29000 in first axis)

    # compute impact parameters, p
    rmax = r_outer[-1]
    ps = calculate_p_values_jax(rmax, N)

    # compute interaction points, z
    zs, shellids, size_zs = populate_z_jax(
        ps, r_inner, r_outer, time_explosion, size_shell
    )

    # init Inup with the photopshere
    Inup = init_Inup_jax(inu, ps, zs[:, 0], iT, r_inner[0])

    # calculate the final intensity
    Inup_i = calc_Inup_jax(
        Inup[:, 1:],
        inu,
        ps[1:],
        zs[1:],
        size_zs[1:],
        shellids[1:],
        line_list_nu,
        time_explosion,
        electron_density,
        exp_tau,
        att_S_ul.ravel(order="F"),
        Jred_lu.ravel(order="F"),
        Jblue_lu.ravel(order="F"),
    )

    # compute the luminosity
    L = 8 * jnp.pi * jnp.pi * jnp.trapezoid(Inup_i, dx=rmax/N, axis=1)    

    return L, Inup_i

class JaxFormalIntegrator:

    def __init__(self, geometry, time_explosion, opacity_state, points):
        
        self.geometry = geometry
        self.time_explosion = time_explosion
        self.opacity_state = opacity_state
        self.points = points

    def formal_integral(
        self, 
        iT,
        frequencies, 
        frequencies_size, # remove?
        att_S_ul, 
        Jred_lu, 
        Jblue_lu, 
        tau_sobolevs,
        electron_densities,
        points
    ):
        # get all necessary parameters and put them in jax arrays
        r_inner_i = jnp.array(self.geometry.r_inner)
        r_outer_i = jnp.array(self.geometry.r_outer)
        tau_sobolevs = jnp.array(tau_sobolevs)
        line_list_nu = jnp.array(self.transport_state.opacity_state.line_list_nu)
        frequencies = jnp.array(frequencies)
        att_S_ul = jnp.array(att_S_ul)
        Jred_lu = jnp.array(Jred_lu)
        Jblue_lu = jnp.array(Jblue_lu)
        electron_densities = jnp.array(electron_densities)

        # compute formal integral
        L, Inup = formal_integral_jax(
            r_inner_i,
            r_outer_i,
            self.time_explosion.value,
            tau_sobolevs,
            line_list_nu,
            iT,
            frequencies,
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            electron_densities,
            points
        )

        return L, Inup