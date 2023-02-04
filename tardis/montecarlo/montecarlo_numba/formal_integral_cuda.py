import sys
import numpy as np
from astropy import units as u
from numba import float64, int64, cuda
import math

from tardis.montecarlo.montecarlo_numba.numba_config import SIGMA_THOMSON

C_INV = 3.33564e-11
M_PI = np.arccos(-1)
KB_CGS = 1.3806488e-16
H_CGS = 6.62606957e-27


@cuda.jit
def cuda_formal_integral(
    r_inner,
    r_outer,
    time_explosion,
    line_list_nu,
    iT,
    inu,
    inu_size,
    att_S_ul,
    Jred_lu,
    Jblue_lu,
    tau_sobolev,
    electron_density,
    N,
    L,
    pp,
    exp_tau,
    I_nu,
    z,
    shell_id,
):
    """
    The CUDA version of numba_formal_integral that can run
    on a NVIDIA GPU.

    Parameters
    ----------
    r_inner : array(float64, 1d, C)
        self.geometry.r_inner
    r_outer : array(float64, 1d, C)
        self.geometry.r_outer
    time_explosion: float64
        self.geometry.time_explosion
    line_list_nu : array(float64, 1d, A)
        self.plasma.line_list_nu
    iT : np.float64
    inu : np.float64
    inu_size : int64
    att_S_ul : array(float64, 1d, C)
    Jred_lu : array(float64, 1d, C)
    Jblue_lu : array(float64, 1d, C)
    tau_sobolev : array(float64, 2d, C)
    electron_density : array(float64, 1d, C)
    N : int64
    L : array(float64, 1d, C)
        This is where the results will be stored
    pp : array(float64, 1d, C)
    exp_tau : array(float64, 1d, C)
    I_nu array(floatt64, 2d, C)
    z : array(float64, 2d, C)
    shell_id : array(int64, 2d, C)
    """

    # todo: add all the original todos

    # global read-only values
    size_line, size_shell = tau_sobolev.shape
    size_tau = size_line * size_shell
    R_ph = r_inner[0]  # make sure these are cgs
    R_max = r_outer[size_shell - 1]

    nu_idx = cuda.grid(1)
    # Check to see if CUDA is out of bounds
    if nu_idx >= inu_size:
        return

    # These all have their own version for each thread to avoid race conditions
    I_nu_thread = I_nu[nu_idx]
    z_thread = z[nu_idx]
    shell_id_thread = shell_id[nu_idx]

    offset = 0
    size_z = 0
    idx_nu_start = 0
    direction = 0
    first = 0
    i = 0
    p = 0.0
    nu_start = 0.0
    nu_end = 0.0
    nu = 0.0
    zstart = 0.0
    zend = 0.0
    escat_contrib = 0.0
    escat_op = 0.0
    Jkkp = 0.0
    pexp_tau = 0
    patt_S_ul = 0
    pJred_lu = 0
    pJblue_lu = 0
    pline = 0

    nu = inu[nu_idx]

    # now loop over discrete values along line
    for p_idx in range(1, N):
        escat_contrib = 0.0
        p = pp[p_idx]

        # initialize z intersections for p values
        size_z = populate_z_cuda(
            r_inner, r_outer, time_explosion, p, z_thread, shell_id_thread
        )
        if p <= R_ph:
            I_nu_thread[p_idx] = intensity_black_body_cuda(nu * z_thread[0], iT)
        else:
            I_nu_thread[p_idx] = 0
        nu_start = nu * z_thread[0]
        nu_end = nu * z_thread[1]

        idx_nu_start = line_search_cuda(line_list_nu, nu_start, size_line)
        offset = shell_id_thread[0] * size_line
        # start tracking accumulated e-scattering optical depth
        zstart = time_explosion / C_INV * (1.0 - z_thread[0])
        # Initialize "pointers"
        pline = int(idx_nu_start)
        pexp_tau = int(offset + idx_nu_start)
        patt_S_ul = int(offset + idx_nu_start)
        pJred_lu = int(offset + idx_nu_start)
        pJblue_lu = int(offset + idx_nu_start)

        # flag for first contribution to integration on current p-ray
        first = 1

        # loop over all interactions
        for i in range(size_z - 1):
            escat_op = electron_density[int(shell_id_thread[i])] * SIGMA_THOMSON
            nu_end = (
                nu * z_thread[i + 1]
            )  # +1 is the offset as the original is from z[1:]

            nu_end_idx = line_search_cuda(
                line_list_nu, nu_end, len(line_list_nu)
            )

            for _ in range(max(nu_end_idx - pline, 0)):

                # calculate e-scattering optical depth to next resonance point
                zend = time_explosion / C_INV * (1.0 - line_list_nu[pline] / nu)
                if first == 1:
                    # first contribution to integration
                    # NOTE: this treatment of I_nu_b (given
                    #   by boundary conditions) is not in Lucy 1999;
                    #   should be re-examined carefully
                    escat_contrib += (
                        (zend - zstart)
                        * escat_op
                        * (Jblue_lu[pJblue_lu] - I_nu_thread[p_idx])
                    )
                    first = 0
                else:
                    # Account for e-scattering, c.f. Eqs 27, 28 in Lucy 1999
                    Jkkp = 0.5 * (Jred_lu[pJred_lu] + Jblue_lu[pJblue_lu])
                    escat_contrib += (
                        (zend - zstart) * escat_op * (Jkkp - I_nu_thread[p_idx])
                    )
                    # this introduces the necessary offset of one element between
                    # pJblue_lu and pJred_lu
                    pJred_lu += 1
                I_nu_thread[p_idx] += escat_contrib
                # // Lucy 1999, Eq 26
                I_nu_thread[p_idx] *= exp_tau[pexp_tau]
                I_nu_thread[p_idx] += att_S_ul[patt_S_ul]

                # // reset e-scattering opacity
                escat_contrib = 0.0
                zstart = zend

                pline += 1
                pexp_tau += 1
                patt_S_ul += 1
                pJblue_lu += 1

            # calculate e-scattering optical depth to grid cell boundary

            Jkkp = 0.5 * (Jred_lu[pJred_lu] + Jblue_lu[pJblue_lu])
            zend = time_explosion / C_INV * (1.0 - nu_end / nu)
            escat_contrib += (
                (zend - zstart) * escat_op * (Jkkp - I_nu_thread[p_idx])
            )
            zstart = zend

            # advance pointers
            direction = int(
                (shell_id_thread[i + 1] - shell_id_thread[i]) * size_line
            )
            pexp_tau += direction
            patt_S_ul += direction
            pJred_lu += direction
            pJblue_lu += direction

        I_nu_thread[p_idx] *= p
    L[nu_idx] = (
        8 * M_PI * M_PI * trapezoid_integration_cuda(I_nu_thread, R_max / N)
    )


class CudaFormalIntegrator(object):
    """
    Helper class for performing the formal integral
    with CUDA.
    """

    def __init__(self, geometry, model, plasma, points=1000):
        self.geometry = geometry
        self.model = model
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
        Simple wrapper for the CUDA implementation of the formal integral
        """
        L = np.zeros(inu_size, dtype=np.float64)  # array(float64, 1d, C)
        # global read-only values
        size_line, size_shell = tau_sobolev.shape  # int64, int64
        size_tau = size_line * size_shell  # int64

        pp = np.zeros(N, dtype=np.float64)  # array(float64, 1d, C)
        exp_tau = np.zeros(size_tau, dtype=np.float64)  # array(float64, 1d, C)
        exp_tau = np.exp(-tau_sobolev.T.ravel())  # array(float64, 1d, C)
        pp[::] = calculate_p_values(
            self.geometry.r_outer[size_shell - 1], N
        )  # array(float64, 1d, C)

        I_nu = np.zeros(
            (inu_size, N), dtype=np.float64
        )  # array(float64, 1d, C)
        z = np.zeros(
            (inu_size, 2 * size_shell), dtype=np.float64
        )  # array(float64, 1d, C)
        shell_id = np.zeros(
            (inu_size, 2 * size_shell), dtype=np.int64
        )  # array(int64, 1d, C)
        THREADS_PER_BLOCK = 32
        blocks_per_grid = (inu_size // THREADS_PER_BLOCK) + 1

        cuda_formal_integral[blocks_per_grid, THREADS_PER_BLOCK](
            self.geometry.r_inner,
            self.geometry.r_outer,
            self.model.time_explosion,
            self.plasma.line_list_nu,
            iT.value,
            inu.value,
            inu_size,
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            tau_sobolev,
            electron_density,
            N,
            L,
            pp,
            exp_tau,
            I_nu,
            z,
            shell_id,
        )

        return L, I_nu


@cuda.jit(device=True)
def populate_z_cuda(r_inner, r_outer, time_explosion, p, oz, oshell_id):
    """
    Calculate p line intersections

    This function calculates the intersection points of the p-line with
    each shell

    Parameters
    ----------
    r_inner : array(float64, 1d, C)
    r_outer : array(float64, 1d, C)
    p : float64
        distance of the integration line to the center
    oz : array(float64, 1d, C)
        will be set with z values. the array is truncated
        by the value `1`.
    oshell_id : array(int64, 1d, C)
        will be set with the corresponding shell_ids

    Returns
    -------
    int64
    """
    N = len(r_inner)
    inv_t = 1 / time_explosion
    z = 0
    offset = N

    if p <= r_inner[0]:
        # intersect the photosphere
        for i in range(N):
            oz[i] = 1 - calculate_z_cuda(r_outer[i], p, inv_t)
            oshell_id[i] = i
        return N
    else:
        # no intersection with photosphere
        # that means we intersect each shell twice
        for i in range(N):
            z = calculate_z_cuda(r_outer[i], p, inv_t)
            if z == 0:
                continue
            if offset == N:
                offset = i
            # calculate the index in the resulting array
            i_low = N - i - 1  # the far intersection with the shell
            i_up = N + i - 2 * offset  # the nearer intersection with the shell

            # setting the arrays
            oz[i_low] = 1 + z
            oshell_id[i_low] = i
            oz[i_up] = 1 - z
            oshell_id[i_up] = i
        return 2 * (N - offset)


@cuda.jit(device=True)
def calculate_z_cuda(r, p, inv_t):
    """
    Calculate distance to p line

    Calculate half of the length of the p-line inside a shell
    of radius r in terms of unit length (c * t_exp).
    If shell and p-line do not intersect, return 0.

    Parameters
    ----------
    r : float64
        radius of the shell
    p : float64
        distance of the p-line to the center of the supernova
    inv_t : float64
        inverse time_explosio is needed to norm to unit-length

    Returns
    -------
    float64
    """
    if r > p:
        return math.sqrt(r * r - p * p) * C_INV * inv_t
    else:
        return 0.0


class BoundsError(IndexError):
    """
    Used to check bounds in reverse
    binary search
    """

    pass


@cuda.jit(device=True)
def line_search_cuda(nu, nu_insert, number_of_lines):
    """
    Insert a value in to an array of line frequencies

    Parameters
    ----------
    nu : (array) line frequencies
    nu_insert : (int) value of nu key
    number_of_lines : (int) number of lines in the line list

    Returns
    -------
    int
        index of the next line to the red.
        If the key value is redder than
        the reddest line returns number_of_lines.
    """
    imin = 0
    imax = number_of_lines - 1
    if nu_insert > nu[imin]:
        result = imin
    elif nu_insert < nu[imax]:
        result = imax + 1
    else:
        result = reverse_binary_search_cuda(nu, nu_insert, imin, imax)
        result = result + 1
    return result


# Credit for this computation is https://github.com/numba/numba/blob/3fd158f79a12ac5276bc5a72c2404464487c91f0/numba/np/arraymath.py#L3542
@cuda.jit(device=True)
def reverse_binary_search_cuda(x, x_insert, imin, imax):
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
    x : np.ndarray(np.float64, 1d, C)
    x_insert : float64
    imin : int
        Lower bound
    imax : int
        Upper bound

    Returns
    -------
    np.int64
        Location of insertion
    """
    if (x_insert > x[imin]) or (x_insert < x[imax]):
        raise BoundsError  # check
    arr = x[::-1]
    n = len(arr)
    lo = 0
    hi = n
    while hi > lo:
        mid = (lo + hi) >> 1
        if arr[mid] <= x_insert:
            # mid is too low of an index, go higher
            lo = mid + 1
        else:
            # mid is too high of an index, go down some
            hi = mid

    return len(x) - 1 - lo


@cuda.jit(device=True)
def trapezoid_integration_cuda(arr, dx):
    """
    Computes the approximation of the
    trapezoidal integration of the array.

    Parameters
    ----------
    arr : (array(float64, 1d, C)
    dx : np.float64
    """

    result = arr[0] + arr[-1]

    for x in range(1, len(arr) - 1):
        result += arr[x] * 2.0

    return result * (dx / 2.0)


@cuda.jit(device=True)
def intensity_black_body_cuda(nu, temperature):
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
    return coefficient * nu * nu * nu / (math.exp(H_CGS * nu * beta_rad) - 1)


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
