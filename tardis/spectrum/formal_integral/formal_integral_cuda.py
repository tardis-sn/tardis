import math
import numpy as np
from numba import cuda

from tardis.spectrum.formal_integral.base import (
    C_INV,
    KB_CGS,
    H_CGS,
    calculate_impact_parameters,
)
from tardis.transport.montecarlo.configuration.constants import SIGMA_THOMSON


class BoundsError(IndexError):
    """
    Used to check bounds in reverse
    binary search
    """


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


@cuda.jit
def cuda_vector_integrator(luminosity_density, intensities_nu, N, radius_max):
    """
    The CUDA Vectorized integrator over second axis

    Parameters
    ----------
    luminosity_density : array(float64, 1d, C)
        Output Array
    intensities_nu : array(float64, 2d, C)
        Input Array
    N : int64
    radius_max : float64

    """
    frequencies_idx = cuda.grid(1)
    luminosity_density[frequencies_idx] = (
        8
        * np.pi
        * np.pi
        * trapezoid_integration_cuda(
            intensities_nu[frequencies_idx], radius_max / N
        )
    )


@cuda.jit(device=True)
def calculate_intersection_point_cuda(
    radius, impact_parameter, inverse_time_explosion
):
    """
    Calculate distance to p line

    Calculate half of the length of the p-line inside a shell
    of radius radius in terms of unit length (c * t_exp).
    If shell and p-line do not intersect, return 0.

    Parameters
    ----------
    radius : float64
        radius of the shell
    impact_parameter : float64
        distance of the p-line to the center of the supernova
    inverse_time_explosion : float64
        inverse time_explosion is needed to norm to unit-length

    Returns
    -------
    float64
    """
    if radius > impact_parameter:
        return (
            math.sqrt(radius * radius - impact_parameter * impact_parameter)
            * C_INV
            * inverse_time_explosion
        )
    else:
        return 0.0


@cuda.jit(device=True)
def populate_intersection_points_cuda(
    radii_inner, radii_outer, time_explosion, impact_parameter, oz, oshell_id
):
    """
    Calculate p line intersections

    This function calculates the intersection points of the p-line with
    each shell

    Parameters
    ----------
    radii_inner : array(float64, 1d, C)
    radii_outer : array(float64, 1d, C)
    time_explosion : float64
    impact_parameter : float64
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
    N = len(radii_inner)
    inverse_time_explosion = 1 / time_explosion
    intersection_point = 0
    offset = N

    if impact_parameter <= radii_inner[0]:
        # intersect the photosphere
        for i in range(N):
            oz[i] = 1 - calculate_intersection_point_cuda(
                radii_outer[i], impact_parameter, inverse_time_explosion
            )
            oshell_id[i] = i
        return N
    else:
        # no intersection with photosphere
        # that means we intersect each shell twice
        for i in range(N):
            intersection_point = calculate_intersection_point_cuda(
                radii_outer[i], impact_parameter, inverse_time_explosion
            )
            if intersection_point == 0:
                continue
            if offset == N:
                offset = i
            # calculate the index in the resulting array
            i_low = N - i - 1  # the far intersection with the shell
            i_up = N + i - 2 * offset  # the nearer intersection with the shell

            # setting the arrays
            oz[i_low] = 1 + intersection_point
            oshell_id[i_low] = i
            oz[i_up] = 1 - intersection_point
            oshell_id[i_up] = i
        return 2 * (N - offset)


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
def line_search_cuda(frequencies, frequencies_insert, number_of_lines):
    """
    Insert a value in to an array of line frequencies

    Parameters
    ----------
    frequencies : (array) line frequencies
    frequencies_insert : (int) value of frequencies key
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
    if frequencies_insert > frequencies[imin]:
        result = imin
    elif frequencies_insert < frequencies[imax]:
        result = imax + 1
    else:
        result = reverse_binary_search_cuda(
            frequencies, frequencies_insert, imin, imax
        )
        result = result + 1
    return result


calculate_impact_parameters = cuda.jit(calculate_impact_parameters, device=True)


@cuda.jit(device=True)
def intensity_black_body_cuda(frequency, temperature):
    """
    Calculate the blackbody intensity.

    Parameters
    ----------
    frequency : float64
        frequency
    temperature : float64
        Temperature

    Returns
    -------
    float64
    """
    if frequency == 0:
        return np.nan  # to avoid ZeroDivisionError
    beta_rad = 1 / (KB_CGS * temperature)
    coefficient = 2 * H_CGS * C_INV * C_INV
    return (
        coefficient
        * frequency
        * frequency
        * frequency
        / (math.exp(H_CGS * frequency * beta_rad) - 1)
    )


@cuda.jit
def cuda_formal_integral(
    radii_inner,
    radii_outer,
    time_explosion,
    line_list_frequencies,
    iT,
    interpolated_frequencies,
    interpolated_frequencies_size,
    att_S_ul,
    Jred_lu,
    Jblue_lu,
    tau_sobolev,
    electron_density,
    N,
    impact_parameter_array,
    exp_tau,
    intensities_nu,
    intersection_point,
    shell_id,
):
    """
    The CUDA version of numba_formal_integral that can run
    on a NVIDIA GPU.

    Parameters
    ----------
    radii_inner : array(float64, 1d, C)
        inner radius of each shell
    radii_outer : array(float64, 1d, C)
        outer radius of each shell
    time_explosion: float64
        geometrical explosion time
    line_list_frequencies : array(float64, 1d, A)
        List of line transition frequencies
    iT : np.float64
        interpolated temperature in cgs units
    interpolated_frequencies : np.float64
        interpolated frequencies in cgs units
    interpolated_frequencies_size : int64
        size of interpolated_frequencies array
    att_S_ul : array(float64, 1d, C)
        attenuated source function
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
    impact_parameter_array : array(float64, 1d, C)
        Impact parameter arrays
    exp_tau : array(float64, 1d, C)
        $\\exp{-tau}$ array to speed up computation
    intensities_nu : array(float64, 2d, C)
        Radiative intensity per unit frequency per impact parameter
    intersection_point : array(float64, 2d, C)
        Ray intersections with the shells
    shell_id : array(int64, 2d, C)
        List of shells for each thread
    """
    # global read-only values
    size_line, size_shell = tau_sobolev.shape
    R_ph = radii_inner[0]  # make sure these are cgs

    frequencies_idx, impact_parameter_idx = cuda.grid(2)  # 2D Cuda Grid, nu x p
    impact_parameter_idx += 1  # Because the iteration starts at 1
    # Check to see if CUDA is out of bounds
    if frequencies_idx >= interpolated_frequencies_size:
        return

    if impact_parameter_idx >= N:
        return

    # These all have their own version for each thread to avoid race conditions
    intensities_nu_thread = intensities_nu[frequencies_idx]
    intersection_point_thread = intersection_point[impact_parameter_idx]
    shell_id_thread = shell_id[impact_parameter_idx]

    offset = 0
    size_z = 0
    idx_frequencies_start = 0
    direction = 0
    first = 0
    i = 0
    impact_parameter = 0.0
    frequency_start = 0.0
    frequency_end = 0.0
    frequency = 0.0
    intersection_point_start = 0.0
    intersection_point_end = 0.0
    escat_contrib = 0.0
    escat_op = 0.0
    Jkkp = 0.0
    pexp_tau = 0
    patt_S_ul = 0
    pJred_lu = 0
    pJblue_lu = 0
    pline = 0

    frequency = interpolated_frequencies[frequencies_idx]
    escat_contrib = 0.0
    impact_parameter = impact_parameter_array[impact_parameter_idx]

    # initialize z intersections for p values
    size_z = populate_intersection_points_cuda(
        radii_inner,
        radii_outer,
        time_explosion,
        impact_parameter,
        intersection_point_thread,
        shell_id_thread,
    )
    if impact_parameter <= R_ph:
        intensities_nu_thread[impact_parameter_idx] = intensity_black_body_cuda(
            frequency * intersection_point_thread[0], iT
        )
    else:
        intensities_nu_thread[impact_parameter_idx] = 0
    frequency_start = frequency * intersection_point_thread[0]
    frequency_end = frequency * intersection_point_thread[1]

    idx_frequencies_start = line_search_cuda(
        line_list_frequencies, frequency_start, size_line
    )
    offset = shell_id_thread[0] * size_line
    # start tracking accumulated e-scattering optical depth
    intersection_point_start = (
        time_explosion / C_INV * (1.0 - intersection_point_thread[0])
    )
    # Initialize "pointers"
    pline = int(idx_frequencies_start)
    pexp_tau = int(offset + idx_frequencies_start)
    patt_S_ul = int(offset + idx_frequencies_start)
    pJred_lu = int(offset + idx_frequencies_start)
    pJblue_lu = int(offset + idx_frequencies_start)

    first = 1

    # loop over all interactions
    for i in range(size_z - 1):
        escat_op = electron_density[int(shell_id_thread[i])] * SIGMA_THOMSON
        frequency_end = (
            frequency * intersection_point_thread[i + 1]
        )  # +1 is the offset as the original is from z[1:]

        frequency_end_idx = line_search_cuda(
            line_list_frequencies, frequency_end, len(line_list_frequencies)
        )
        for _ in range(max(frequency_end_idx - pline, 0)):
            # calculate e-scattering optical depth to next resonance point
            intersection_point_end = (
                time_explosion
                / C_INV
                * (1.0 - line_list_frequencies[pline] / frequency)
            )  # check

            if first == 1:
                escat_contrib += (
                    (intersection_point_end - intersection_point_start)
                    * escat_op
                    * (
                        Jblue_lu[pJblue_lu]
                        - intensities_nu_thread[impact_parameter_idx]
                    )
                )
                first = 0
            else:
                # Account for e-scattering, c.f. Eqs 27, 28 in Lucy 1999
                Jkkp = 0.5 * (Jred_lu[pJred_lu] + Jblue_lu[pJblue_lu])
                escat_contrib += (
                    (intersection_point_end - intersection_point_start)
                    * escat_op
                    * (Jkkp - intensities_nu_thread[impact_parameter_idx])
                )
                # this introduces the necessary ffset of one element between
                # pJblue_lu and pJred_lu
                pJred_lu += 1
            intensities_nu_thread[impact_parameter_idx] += escat_contrib
            # // Lucy 1999, Eq 26
            intensities_nu_thread[impact_parameter_idx] *= exp_tau[pexp_tau]
            intensities_nu_thread[impact_parameter_idx] += att_S_ul[patt_S_ul]

            # // reset e-scattering opacity
            escat_contrib = 0
            intersection_point_start = intersection_point_end

            pline += 1
            pexp_tau += 1
            patt_S_ul += 1
            pJblue_lu += 1

        # calculate e-scattering optical depth to grid cell boundary

        Jkkp = 0.5 * (Jred_lu[pJred_lu] + Jblue_lu[pJblue_lu])
        intersection_point_end = (
            time_explosion / C_INV * (1.0 - frequency_end / frequency)
        )
        escat_contrib += (
            (intersection_point_end - intersection_point_start)
            * escat_op
            * (Jkkp - intensities_nu_thread[impact_parameter_idx])
        )
        intersection_point_start = intersection_point_end

        # advance pointers
        direction = int(
            (shell_id_thread[i + 1] - shell_id_thread[i]) * size_line
        )
        pexp_tau += direction
        patt_S_ul += direction
        pJred_lu += direction
        pJblue_lu += direction

    intensities_nu_thread[impact_parameter_idx] *= impact_parameter


def calculate_impact_parameter_values(radius_max, N):
    """
    Calculates the impact parameter values of N

    Parameters
    ----------
    radius_max : float64
    N : int64

    Returns
    -------
    float64
    """
    return np.arange(N).astype(np.float64) * radius_max / (N - 1)


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
        interpolated_frequencies,
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
        interpolated_frequencies_size = len(interpolated_frequencies)

        impact_parameters_array = np.zeros(
            N, dtype=np.float64
        )  # array(float64, 1d, C)
        exp_tau = np.zeros(size_tau, dtype=np.float64)  # array(float64, 1d, C)
        exp_tau = np.exp(-tau_sobolev.T.ravel())  # array(float64, 1d, C)
        impact_parameters_array[::] = calculate_impact_parameter_values(
            self.geometry.radii_outer[size_shell - 1], N
        )  # array(float64, 1d, C)
        intersection_point = np.zeros(
            (N, 2 * size_shell), dtype=np.float64
        )  # array(float64, 2d, C)
        shell_id = np.zeros(
            (N, 2 * size_shell), dtype=np.int64
        )  # array(int64, 2d, C)

        # These get separate names since they'll be copied back
        # These are device objects stored on the GPU
        # for the Luminosity density and Radiative intensity
        d_lum_density = cuda.device_array(
            (interpolated_frequencies_size,), dtype=np.float64
        )
        d_intensities_nu = cuda.device_array(
            (interpolated_frequencies_size, N), dtype=np.float64
        )

        # Copy these arrays to the device, we don't need them again
        # But they must be initialized with zeros
        intersection_point = cuda.to_device(intersection_point)
        shell_id = cuda.to_device(shell_id)
        impact_parameters_array = cuda.to_device(impact_parameters_array)
        exp_tau = cuda.to_device(exp_tau)
        radii_inner = cuda.to_device(self.geometry.radii_inner)
        radii_outer = cuda.to_device(self.geometry.radii_outer)
        line_list_frequencies = cuda.to_device(
            self.plasma.line_list_frequencies
        )
        interpolated_frequencies = cuda.to_device(
            interpolated_frequencies.value
        )
        att_S_ul = cuda.to_device(att_S_ul)
        Jred_lu = cuda.to_device(Jred_lu)
        Jblue_lu = cuda.to_device(Jblue_lu)
        tau_sobolev = cuda.to_device(tau_sobolev)
        electron_density = cuda.to_device(electron_density)

        # Thread/Block Allocation, this seems to work
        THREADS_PER_BLOCK_NU = 32
        THREADS_PER_BLOCK_P = 16
        blocks_per_grid_nu = (
            interpolated_frequencies_size // THREADS_PER_BLOCK_NU
        ) + 1
        blocks_per_grid_p = ((N - 1) // THREADS_PER_BLOCK_P) + 1

        cuda_formal_integral[
            (blocks_per_grid_nu, blocks_per_grid_p),
            (THREADS_PER_BLOCK_NU, THREADS_PER_BLOCK_P),
        ](
            radii_inner,
            radii_outer,
            self.time_explosion,
            line_list_frequencies,
            iT.value,
            interpolated_frequencies,
            interpolated_frequencies_size,
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            tau_sobolev,
            electron_density,
            N,
            impact_parameters_array,
            exp_tau,
            d_intensities_nu,
            intersection_point,
            shell_id,
        )

        radius_max = self.geometry.radii_outer[size_shell - 1]
        cuda_vector_integrator[blocks_per_grid_nu, THREADS_PER_BLOCK_NU](
            d_lum_density, d_intensities_nu, N, radius_max
        )
        lum_density = d_lum_density.copy_to_host()
        intensities_nu = d_intensities_nu.copy_to_host()

        return lum_density, intensities_nu
