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

@njit(**njit_dict_no_parallel)
def increment_escat_contrib(
        escat_contrib,
        first,
        pJred_lu,
        zend,
        zstart,
        escat_op,
        Jblue_lu, 
        Jred_lu,
        I_nu_p
    ):

    if first == 1:
        # first contribution to integration
        # NOTE: this treatment of I_nu_b (given
        #   by boundary conditions) is not in Lucy 1999;
        #   should be re-examined carefully
        escat_contrib += (
            (zend - zstart)
            * escat_op
            * (Jblue_lu - I_nu_p)
        )
        first = 0
    else:
        # Account for e-scattering, c.f. Eqs 27, 28 in Lucy 1999
        Jkkp = 0.5 * (Jred_lu + Jblue_lu)
        escat_contrib += (
            (zend - zstart) * escat_op * (Jkkp - I_nu_p)
        )
        # this introduces the necessary ffset of one element between
        # pJblue_lu and pJred_lu
        pJred_lu += 1

    return escat_contrib, first, pJred_lu

# initialize I_nu_p with the blackbody intensities (as necessary)
@njit(**njit_dict)
def init_Inup_numba(inu, ps, iT, Rph, N, size_shell, geometry, time_explosion):
    inu_size = len(inu)

    I_nu_p = np.zeros((inu_size, N), dtype=np.float64)
    zs = np.zeros((N, 2 * size_shell), dtype=np.float64)
    shell_ids = np.zeros((N, 2 * size_shell), dtype=np.int64)
    size_zs = np.zeros(N, dtype=np.int64)

    for nu_idx in prange(inu_size):
        I_nu = I_nu_p[nu_idx]
        nu = inu[nu_idx]
        for p_idx in range(1, N):
            p = ps[p_idx]

            size_z = populate_z(geometry, time_explosion, p, zs[p_idx], shell_ids[p_idx])
            size_zs[p_idx] = size_z

            if p <= Rph:
                I_nu[p_idx] = intensity_black_body(nu * zs[p_idx][0], iT)
            else:
                I_nu[p_idx] = 0
    
    return I_nu_p, zs, shell_ids, size_zs


@njit(**njit_dict_no_parallel)
def numba_formal_integral(
    geometry,
    time_explosion,
    plasma,
    iT,
    inu,
    inu_size,
    att_S_ul,
    Jred_lu,
    Jblue_lu,
    tau_sobolev,
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
    # TODO: add all the original todos
    # Initialize the output which is shared among threads
    L = np.zeros(inu_size, dtype=np.float64)
    # global read-only values
    size_line, size_shell = tau_sobolev.shape
    size_tau = size_line * size_shell
    R_ph = geometry.r_inner[0]  # make sure these are cgs
    R_max = geometry.r_outer[size_shell - 1]
    pp = np.zeros(N, dtype=np.float64)  # check
    exp_tau = np.zeros(size_tau, dtype=np.float64)
    exp_tau = np.exp(-tau_sobolev.T.ravel())  # maybe make this 2D?
    pp[::] = calculate_p_values(R_max, N)
    line_list_nu = plasma.line_list_nu
    # done with instantiation
    # now loop over wavelength in spectrum
    # I_nu_p = np.zeros((inu_size, N), dtype=np.float64)

    I_nu_p, zs, shell_ids, size_zs = init_Inup_numba(
        inu, pp, iT, R_ph, N, size_shell, geometry, time_explosion
    )

    for nu_idx in prange(inu_size):
        I_nu = I_nu_p[nu_idx]
        first = 0

        nu = inu[nu_idx]
        # now loop over discrete values along line
        for p_idx in range(1, N):
            escat_contrib = 0
            p = pp[p_idx]
            z = zs[p_idx]
            shell_id = shell_ids[p_idx]
            size_z = size_zs[p_idx]

            # init black_body moved up

            # find first contributing lines
            nu_start = nu * z[0]
            nu_end = nu * z[1]
            idx_nu_start = line_search(plasma.line_list_nu, nu_start, size_line)
            offset = shell_id[0] * size_line
            # start tracking accumulated e-scattering optical depth
            zstart = time_explosion / C_INV * (1.0 - z[0])
            # Initialize "pointers"
            pline = int(idx_nu_start)
            pexp_tau = int(offset + idx_nu_start)
            patt_S_ul = int(offset + idx_nu_start)
            pJred_lu = int(offset + idx_nu_start)
            pJblue_lu = int(offset + idx_nu_start)

            # flag for first contribution to integration on current p-ray
            first = 1
            nu_ends = nu * z[1:]
            nu_ends_idxs = size_line - np.searchsorted(
                line_list_nu[::-1], nu_ends, side="right"
            )
            # loop over all interactions
            for i in range(size_z - 1):
                escat_op = electron_density[int(shell_id[i])] * SIGMA_THOMSON
                nu_end = nu_ends[i]
                nu_end_idx = nu_ends_idxs[i]
                for _ in range(max(nu_end_idx - pline, 0)):
                    # calculate e-scattering optical depth to next resonance point
                    zend = (
                        time_explosion
                        / C_INV
                        * (1.0 - line_list_nu[pline] / nu)
                    )  # check

                    Jblue_lu_i = Jblue_lu[pJblue_lu]
                    Jred_lu_i = Jred_lu[pJred_lu]
                    I_nu_i = I_nu[p_idx]
                    escat_contrib, first, pJred_lu = increment_escat_contrib(
                        escat_contrib, 
                        first,
                        pJred_lu,
                        zend,
                        zstart,
                        escat_op,
                        Jblue_lu_i,
                        Jred_lu_i,
                        I_nu_i
                    )

                    I_nu[p_idx] += escat_contrib
                    # // Lucy 1999, Eq 26
                    I_nu[p_idx] *= exp_tau[pexp_tau]
                    I_nu[p_idx] += att_S_ul[patt_S_ul]

                    # // reset e-scattering opacity
                    escat_contrib = 0
                    zstart = zend

                    pline += 1
                    pexp_tau += 1
                    patt_S_ul += 1
                    pJblue_lu += 1

                # calculate e-scattering optical depth to grid cell boundary

                Jkkp = 0.5 * (Jred_lu[pJred_lu] + Jblue_lu[pJblue_lu])
                zend = time_explosion / C_INV * (1.0 - nu_end / nu)  # check
                escat_contrib += (
                    (zend - zstart) * escat_op * (Jkkp - I_nu[p_idx])
                )
                zstart = zend

                # advance pointers
                direction = int((shell_id[i + 1] - shell_id[i]) * size_line)
                pexp_tau += direction
                patt_S_ul += direction
                pJred_lu += direction
                pJblue_lu += direction
            I_nu[p_idx] *= p
        L[nu_idx] = 8 * np.pi * np.pi * trapezoid_integration(I_nu, R_max / N)

    return L, I_nu_p


# integrator_spec = [
#    ("model", NumbaModel.class_type.instance_type),
#    ("plasma", OpacityState.class_type.instance_type),
#    ("points", int64),
# ]


# @jitclass(integrator_spec)
class NumbaFormalIntegrator:
    """
    Helper class for performing the formal integral
    with numba.
    """

    def __init__(self, geometry, time_explosion, plasma, points=1000):
        self.geometry = geometry
        self.time_explosion = time_explosion
        self.plasma = plasma
        self.points = points

    def formal_integral(
        self,
        iT,
        inu,
        inu_size,
        att_S_ul,
        Jred_lu,
        Jblue_lu,
        tau_sobolev,
        electron_density,
        N,
    ):
        """
        Simple wrapper for the numba implementation of the formal integral
        """
        return numba_formal_integral(
            self.geometry,
            self.time_explosion,
            self.plasma,
            iT,
            inu,
            inu_size,
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            tau_sobolev,
            electron_density,
            N,
        )
