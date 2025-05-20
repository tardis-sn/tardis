import math

import numpy as np
from numba import cuda

from tardis.transport.montecarlo.configuration.constants import SIGMA_THOMSON

C_INV = 3.33564e-11
PI = np.pi
KB_CGS = 1.3806488e-16
H_CGS = 6.62606957e-27


@cuda.jit
def cuda_vector_integrator(L, I_nu, N, R_max):
    """
    The CUDA Vectorized integrator over second axis

    Parameters
    ----------
    L : array(float64, 1d, C)
        Output Array
    I_nu : array(floagt64, 2d, C)
        Input Array
    N : int64
    R_max : float64

    """
    nu_idx = cuda.grid(1)
    L[nu_idx] = (
        8 * PI * PI * trapezoid_integration_cuda(I_nu[nu_idx], R_max / N)
    )


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
        inner radius of each shell
    r_outer : array(float64, 1d, C)
        outer radius of each shell
    time_explosion: float64
        geometrical explosion time
    line_list_nu : array(float64, 1d, A)
        List of line transition frequencies
    iT : np.float64
        interpolated temperture in cgs units
    inu : np.float64
        interpolated frequencies in cgs units
    inu_size : int64
        size of inu array
    att_S_ul : array(float64, 1d, C)
        attentuated source function
    Jred_lu : array(float64, 1d, C)
        J estimator from red end of the line from lower to upper level
    Jblue_lu : array(float64, 1d, C)
        J estimator from blue end of the line from lower to upper level
    tau_sobolev : array(float64, 2d, C)
        Sobolev Optical depth for each line in each shell
    electron_density : array(float64, 1d, C)
        electron density in each shell
    N : int64
        Number of impact parameter values (p)
    L : array(float64, 1d, C)
        Luminosity density at each frequency
        This is where the results will be stored
    pp : array(float64, 1d, C)
        Impact parameter arrays
    exp_tau : array(float64, 1d, C)
        $\\exp{-tau}$ array to speed up computation
    I_nu array(floatt64, 2d, C)
        Radiative intensity per unit frequency per impact parameter
    z : array(float64, 2d, C)
        Ray intersections with the shells
    shell_id : array(int64, 2d, C)
        List of shells for each thread
    """
    # global read-only values
    size_line, size_shell = tau_sobolev.shape
    R_ph = r_inner[0]  # make sure these are cgs

    nu_idx, p_idx = cuda.grid(2)  # 2D Cuda Grid, nu x p
    p_idx += 1  # Because the iteration starts at 1
    # Check to see if CUDA is out of bounds
    if nu_idx >= inu_size:
        return

    if p_idx >= N:
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

    first = 1

    # loop over all interactions
    for i in range(size_z - 1):
        escat_op = electron_density[int(shell_id_thread[i])] * SIGMA_THOMSON
        nu_end = (
            nu * z_thread[i + 1]
        )  # +1 is the offset as the original is from z[1:]

        nu_end_idx = line_search_cuda(line_list_nu, nu_end, len(line_list_nu))
        for _ in range(max(nu_end_idx - pline, 0)):
            # calculate e-scattering optical depth to next resonance point
            zend = (
                time_explosion / C_INV * (1.0 - line_list_nu[pline] / nu)
            )  # check

            if first == 1:
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
                # this introduces the necessary ffset of one element between
                # pJblue_lu and pJred_lu
                pJred_lu += 1
            I_nu_thread[p_idx] += escat_contrib
            # // Lucy 1999, Eq 26
            I_nu_thread[p_idx] *= exp_tau[pexp_tau]
            I_nu_thread[p_idx] += att_S_ul[patt_S_ul]

            # // reset e-scattering opacity
            escat_contrib = 0
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


class CudaFormalIntegrator:
    """
    Helper class for performing the formal integral
    with CUDA.
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
        Simple wrapper for the CUDA implementation of the formal integral
        """
        # global read-only values
        size_line, size_shell = tau_sobolev.shape  # int64, int64
        size_tau = size_line * size_shell  # int64

        pp = np.zeros(N, dtype=np.float64)  # array(float64, 1d, C)
        exp_tau = np.zeros(size_tau, dtype=np.float64)  # array(float64, 1d, C)
        exp_tau = np.exp(-tau_sobolev.T.ravel())  # array(float64, 1d, C)
        pp[::] = calculate_p_values(
            self.geometry.r_outer[size_shell - 1], N
        )  # array(float64, 1d, C)
        z = np.zeros(
            (inu_size, 2 * size_shell), dtype=np.float64
        )  # array(float64, 1d, C)
        shell_id = np.zeros(
            (inu_size, 2 * size_shell), dtype=np.int64
        )  # array(int64, 1d, C)

        # These get separate names since they'll be copied back
        # These are device objects stored on the GPU
        # for the Luminosity density and Radiative intensity
        d_L = cuda.device_array((inu_size,), dtype=np.float64)
        d_I_nu = cuda.device_array((inu_size, N), dtype=np.float64)

        # Copy these arrays to the device, we don't need them again
        # But they must be initialized with zeros
        z = cuda.to_device(z)
        shell_id = cuda.to_device(shell_id)
        pp = cuda.to_device(pp)
        exp_tau = cuda.to_device(exp_tau)
        r_inner = cuda.to_device(self.geometry.r_inner)
        r_outer = cuda.to_device(self.geometry.r_outer)
        line_list_nu = cuda.to_device(self.plasma.line_list_nu)
        inu = cuda.to_device(inu.value)
        att_S_ul = cuda.to_device(att_S_ul)
        Jred_lu = cuda.to_device(Jred_lu)
        Jblue_lu = cuda.to_device(Jblue_lu)
        tau_sobolev = cuda.to_device(tau_sobolev)
        electron_density = cuda.to_device(electron_density)

        # Thread/Block Allocation, this seems to work
        THREADS_PER_BLOCK_NU = 32
        THREADS_PER_BLOCK_P = 16
        blocks_per_grid_nu = (inu_size // THREADS_PER_BLOCK_NU) + 1
        blocks_per_grid_p = ((N - 1) // THREADS_PER_BLOCK_P) + 1

        cuda_formal_integral[
            (blocks_per_grid_nu, blocks_per_grid_p),
            (THREADS_PER_BLOCK_NU, THREADS_PER_BLOCK_P),
        ](
            r_inner,
            r_outer,
            self.time_explosion,
            line_list_nu,
            iT.value,
            inu,
            inu_size,
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            tau_sobolev,
            electron_density,
            N,
            pp,
            exp_tau,
            d_I_nu,
            z,
            shell_id,
        )

        R_max = self.geometry.r_outer[size_shell - 1]
        cuda_vector_integrator[blocks_per_grid_nu, THREADS_PER_BLOCK_NU](
            d_L, d_I_nu, N, R_max
        )
        L = d_L.copy_to_host()
        I_nu = d_I_nu.copy_to_host()

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
