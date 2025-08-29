import numpy as np
from numba import njit, prange

from tardis.spectrum.formal_integral.base import (
    C_INV,
    BoundsError,
    calculate_p_values,
    intensity_black_body,
)
from tardis.transport.montecarlo import njit_dict, njit_dict_no_parallel
from tardis.transport.montecarlo.configuration.constants import SIGMA_THOMSON


@njit(**njit_dict_no_parallel)
def calculate_z(r, p, inv_t):
    """
    Calculate distance to p line

    Calculate half of the length of the p-line inside a shell
    of radius r in terms of unit length (c * t_exp).
    If shell and p-line do not intersect, return 0.

    Parameters
    ----------
    r : float
        Radius of the shell.
    p : float
        Distance of the p-line to the center of the supernova.
    inv_t : float
        Inverse of the time_explosion, used to normalize to unit length.

    Returns
    -------
    float
        Half the length of the p-line inside the shell, or 0 if no intersection.
    """
    if r > p:
        return np.sqrt(r * r - p * p) * C_INV * inv_t
    else:
        return 0


@njit(**njit_dict_no_parallel)
def populate_z(geometry, time_explosion, p, zs, shell_ids):
    """
    Calculate the intersection points of the p-line with each shell.

    Parameters
    ----------
    geometry : object
        Geometry object containing shell radii.
    time_explosion : float
        Time since explosion (seconds).
    p : float
        Distance of the integration line to the center.
    zs : ndarray
        Output array to be filled with z values.
    shell_ids : ndarray
        Output array to be filled with the corresponding shell IDs.

    Returns
    -------
    int
        Number of intersection points found.
    """
    # abbreviations
    r = geometry.r_outer
    N = len(geometry.r_inner)
    inv_t = 1 / time_explosion
    z = 0
    offset = N

    if p <= geometry.r_inner[0]:
        # intersect the photosphere
        for i in range(N):
            zs[i] = 1 - calculate_z(r[i], p, inv_t)
            shell_ids[i] = i
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

            zs[i_low] = 1 + z
            shell_ids[i_low] = i
            zs[i_up] = 1 - z
            shell_ids[i_up] = i
        return 2 * (N - offset)


@njit(**njit_dict_no_parallel)
def reverse_binary_search(x, x_insert, imin, imax):
    """
    Find the insertion index for a value in an inversely sorted float array.

    Parameters
    ----------
    x : ndarray
        Inversely (largest to lowest) sorted float array.
    x_insert : float
        Value to insert.
    imin : int
        Lower bound index.
    imax : int
        Upper bound index.

    Returns
    -------
    int
        Index of the next boundary to the left.
    """

    if (x_insert > x[imin]) or (x_insert < x[imax]):
        raise BoundsError
    return len(x) - 1 - np.searchsorted(x[::-1], x_insert, side="right")


@njit(**njit_dict_no_parallel)
def line_search(nu, nu_insert, number_of_lines):
    """
    Find the index to insert a value into an array of line frequencies.

    Parameters
    ----------
    nu : ndarray
        Array of line frequencies.
    nu_insert : float
        Value of the frequency to insert.
    number_of_lines : int
        Number of lines in the line list.

    Returns
    -------
    int
        Index of the next line to the red. If the key value is redder than the reddest line, returns number_of_lines.
    """

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
    """
    Integrate an array using the trapezoidal rule.

    Parameters
    ----------
    array : ndarray
        Array of values to integrate.
    h : float
        Step size.

    Returns
    -------
    float
        The integral of the array using the trapezoidal rule.
    """
    return np.trapezoid(array, dx=h)


# numba jit
calculate_p_values = njit(calculate_p_values, **njit_dict_no_parallel)
intensity_black_body = njit(intensity_black_body, **njit_dict_no_parallel)


@njit(**njit_dict)
def setup_formal_integral_inputs(
    frequencies,
    inner_temperature,
    points,
    geometry,
    time_explosion,
    line_list_nu,
    tau_sobolev,
    electron_densities,
):
    """
    Prepare all arrays and values needed for the loops inside the formal integral.

    Parameters
    ----------
    frequencies : ndarray
        Array of frequency values.
    inner_temperature : float
        Inner boundary temperature.
    points : int
        Number of p-rays.
    geometry : object
        Geometry object containing shell radii.
    time_explosion : float
        Time since explosion (seconds).
    line_list_nu : ndarray
        Array of line frequencies.
    tau_sobolev : ndarray
        Sobolev optical depths for each line and shell.
    electron_densities : ndarray
        Electron densities per shell.

    Returns
    -------
    intensities_nu_p : ndarray
        Intensities at each frequency and impact parameter (p).
    ps : ndarray
        Array of p values.
    size_zs : ndarray
        Number of z intersections for each p.
    lines_idx : ndarray
        Initial line indices for each frequency and p.
    lines_idx_offset : ndarray
        Initial offset line indices for each frequency and p.
    nu_ends : ndarray
        Frequency at each z intersection.
    nu_ends_idxs : ndarray
        Indices of line list for each nu_end.
    zstarts : ndarray
        Starting z values for each p-ray.
    escat_ops : ndarray
        Electron scattering opacities for each frequency, p, and z.
    directions : ndarray
        Direction changes for each frequency, p-ray, and z.
    exp_tau : ndarray
        Exponential of negative Sobolev optical depths (flattened).
    """

    n_frequencies = len(frequencies)
    size_line, size_shell = tau_sobolev.shape
    exp_tau = np.exp(-tau_sobolev.T.ravel())
    R_max = geometry.r_outer[size_shell - 1]
    R_photosphere = geometry.r_inner[0]

    ps = np.zeros(points, dtype=np.float64)
    ps[::] = calculate_p_values(R_max, points)

    # prepare arrays for returned parameters
    intensities_nu_p = np.zeros((n_frequencies, points), dtype=np.float64)
    zs = np.zeros((points, 2 * size_shell), dtype=np.float64)
    shell_ids = np.zeros((points, 2 * size_shell), dtype=np.int64)
    size_zs = np.zeros(points, dtype=np.int64)

    lines_idx = np.zeros((n_frequencies, points), dtype=np.int64)
    lines_idx_offset = np.zeros((n_frequencies, points), dtype=np.int64)

    zstarts = np.zeros(points, dtype=np.float64)

    nu_ends = np.zeros(
        (n_frequencies, points, 2 * size_shell - 1), dtype=np.float64
    )
    nu_ends_idxs = np.zeros(
        (n_frequencies, points, 2 * size_shell - 1), dtype=np.int64
    )

    escat_ops = np.zeros(
        (n_frequencies, points, 2 * size_shell), dtype=np.float64
    )
    directions = np.zeros(
        (n_frequencies, points, 2 * size_shell), dtype=np.int64
    )

    # loop over nu and p
    for nu_idx in prange(n_frequencies):
        intensities_nu = intensities_nu_p[nu_idx]
        nu = frequencies[nu_idx]
        for p_idx in range(1, points):
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
            if p <= R_photosphere:
                intensities_nu[p_idx] = intensity_black_body(
                    nu * z[0], inner_temperature
                )
            else:
                intensities_nu[p_idx] = 0

            # compute quantities for line interactions
            nu_start = nu * z[0]

            idx_nu_start = line_search(line_list_nu, nu_start, size_line)
            lines_idx[nu_idx, p_idx] = int(idx_nu_start)
            lines_idx_offset[nu_idx, p_idx] = int(
                idx_nu_start + (shell_id[0] * size_line)
            )

            zstarts[p_idx] = time_explosion / C_INV * (1.0 - z[0])

            nu_ends[nu_idx, p_idx] = nu * z[1:]
            nu_ends_idxs[nu_idx, p_idx] = size_line - np.searchsorted(
                line_list_nu[::-1], nu_ends[nu_idx, p_idx], side="right"
            )

            for i in range(size_z - 1):
                escat_ops[nu_idx, p_idx, i] = (
                    electron_densities[int(shell_ids[p_idx][i])] * SIGMA_THOMSON
                )
                directions[nu_idx, p_idx, i] = int(
                    (shell_ids[p_idx][i + 1] - shell_ids[p_idx][i]) * size_line
                )

    return (
        intensities_nu_p,
        ps,
        size_zs,
        lines_idx,
        lines_idx_offset,
        nu_ends,
        nu_ends_idxs,
        zstarts,
        escat_ops,
        directions,
        exp_tau,
    )


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
    intensities_nu_p,
):
    """
    Increment the electron scattering contribution for the formal integral.

    Parameters
    ----------
    escat_contrib : float
        Current electron scattering contribution.
    first : int
        Flag indicating if this is the first contribution (1 if first, 0 otherwise).
    pJred_lu : int
        Index for Jred_lu.
    zend : float
        Ending z value for the current segment.
    zstart : float
        Starting z value for the current segment.
    escat_op : float
        Electron scattering opacity.
    Jblue_lu : float
        Jblue_lu value for the current segment.
    Jred_lu : float
        Jred_lu value for the current segment.
    intensities_nu_p : float
        Intensity at the current frequency and p.

    Returns
    -------
    escat_contrib : float
        Updated electron scattering contribution.
    first : int
        Updated first flag.
    pJred_lu : int
        Updated Jred_lu index.
    """
    if first == 1:
        # first contribution to integration
        # NOTE: this treatment of I_nu_b (given
        #   by boundary conditions) is not in Lucy 1999;
        #   should be re-examined carefully
        escat_contrib += (
            (zend - zstart) * escat_op * (Jblue_lu - intensities_nu_p)
        )
        first = 0
    else:
        # Account for e-scattering, c.f. Eqs 27, 28 in Lucy 1999
        Jkkp = 0.5 * (Jred_lu + Jblue_lu)
        escat_contrib += (zend - zstart) * escat_op * (Jkkp - intensities_nu_p)
        # this introduces the necessary offset of one element between
        # pJblue_lu and pJred_lu
        pJred_lu += 1

    return escat_contrib, first, pJred_lu


@njit(**njit_dict)
def numba_formal_integral(
    geometry,
    time_explosion,
    plasma,
    inner_temperature,
    frequencies,
    att_S_ul,
    Jred_lu,
    Jblue_lu,
    tau_sobolev,
    electron_density,
    points,
):
    """
    Compute the formal integral.

    Parameters
    ----------
    geometry : object
        Geometry object containing shell radii.
    time_explosion : float
        Time since explosion (seconds).
    plasma : object
        Plasma object containing line list frequencies.
    inner_temperature : float
        Inner boundary temperature.
    frequencies : ndarray
        Array of frequencies.
    att_S_ul : ndarray
        Attenuated source function for each line and shell.
    Jred_lu : ndarray
        Jred_lu values for each line and shell.
    Jblue_lu : ndarray
        Jblue_lu values for each line and shell.
    tau_sobolev : ndarray
        Sobolev optical depths for each line and shell.
    electron_density : ndarray
        Electron densities per shell.
    points : int
        Number of p-ray points.

    Returns
    -------
    L : ndarray
        Integrated luminosities for each frequency.
    intensities_nu_p : ndarray
        Intensities at each p-ray multiplied by p (frequency x p-ray grid).
    """
    # Initialize the output which is shared among threads
    n_frequencies = len(frequencies)
    L = np.zeros(n_frequencies, dtype=np.float64)

    R_max = geometry.r_outer[-1]
    line_list_nu = plasma.line_list_nu

    # prepare all of the quantities needed for the loops
    (
        intensities_nu_p,
        ps,
        size_zs,
        lines_idx,
        lines_idx_offset,
        nu_ends,
        nu_ends_idxs,
        zstarts,
        escat_ops,
        directions,
        exp_tau,
    ) = setup_formal_integral_inputs(
        frequencies,
        inner_temperature,
        points,
        geometry,
        time_explosion,
        line_list_nu,
        tau_sobolev,
        electron_density,
    )

    # now loop over wavelength in spectrum
    for nu_idx in prange(n_frequencies):
        intensities_nu = intensities_nu_p[nu_idx]
        nu = frequencies[nu_idx]

        # now loop over discrete values along line
        for p_idx in range(1, points):
            escat_contrib = 0
            p = ps[p_idx]
            size_z = size_zs[p_idx]

            # find first contributing lines
            zstart = zstarts[p_idx]

            # Initialize "pointers"
            line_idx = lines_idx[nu_idx, p_idx]
            line_idx_offset = lines_idx_offset[nu_idx, p_idx]
            line_Jred_lu_idx = lines_idx_offset[nu_idx, p_idx]

            # flag for first contribution to integration on current p-ray
            first = 1

            # loop over all interactions
            for i in range(size_z - 1):
                escat_op = escat_ops[nu_idx, p_idx, i]
                nu_end = nu_ends[nu_idx, p_idx, i]
                nu_end_idx = nu_ends_idxs[nu_idx, p_idx, i]
                for _ in range(max(nu_end_idx - line_idx, 0)):
                    # calculate e-scattering optical depth to next resonance point
                    zend = (
                        time_explosion
                        / C_INV
                        * (1.0 - line_list_nu[line_idx] / nu)
                    )

                    escat_contrib, first, line_Jred_lu_idx = (
                        increment_escat_contrib(
                            escat_contrib,
                            first,
                            line_Jred_lu_idx,
                            zend,
                            zstart,
                            escat_op,
                            Jblue_lu[line_idx_offset],
                            Jred_lu[line_Jred_lu_idx],
                            intensities_nu[p_idx],
                        )
                    )

                    intensities_nu[p_idx] += escat_contrib
                    # Lucy 1999, Eq 26
                    intensities_nu[p_idx] *= exp_tau[line_idx_offset]
                    intensities_nu[p_idx] += att_S_ul[line_idx_offset]

                    # reset e-scattering opacity
                    escat_contrib = 0
                    zstart = zend

                    line_idx += 1
                    line_idx_offset += 1

                # calculate e-scattering optical depth to grid cell boundary
                Jkkp = 0.5 * (
                    Jred_lu[line_Jred_lu_idx] + Jblue_lu[line_idx_offset]
                )
                zend = time_explosion / C_INV * (1.0 - nu_end / nu)
                escat_contrib += (
                    (zend - zstart) * escat_op * (Jkkp - intensities_nu[p_idx])
                )
                zstart = zend

                # advance "pointers"
                direction = directions[nu_idx, p_idx, i]
                line_idx_offset += direction
                line_Jred_lu_idx += direction
            intensities_nu[p_idx] *= p
        L[nu_idx] = (
            8
            * np.pi
            * np.pi
            * trapezoid_integration(intensities_nu, R_max / points)
        )

    return L, intensities_nu_p


class NumbaFormalIntegrator:
    """
    Helper class for performing the formal integral with Numba.

    Parameters
    ----------
    geometry : object
        Geometry object containing shell radii.
    time_explosion : float
        Time since explosion (seconds).
    plasma : object
        Plasma object containing line list frequencies.
    points : int, optional
        Number of p-ray points (default is 1000).
    """

    def __init__(self, geometry, time_explosion, plasma, points=1000):
        self.geometry = geometry
        self.time_explosion = time_explosion
        self.plasma = plasma
        self.points = points

    def formal_integral(
        self,
        inner_temperature,
        frequencies,
        att_S_ul,
        Jred_lu,
        Jblue_lu,
        tau_sobolev,
        electron_density,
        points,
    ):
        """
        Wrapper for the Numba implementation of the formal integral.

        Parameters
        ----------
        inner_temperature : float
            Inner boundary temperature.
        frequencies : ndarray
            Array of frequency values.
        att_S_ul : ndarray
            Attenuated source function values.
        Jred_lu : ndarray
            Jred_lu values for each line and shell.
        Jblue_lu : ndarray
            Jblue_lu values for each line and shell.
        tau_sobolev : ndarray
            Sobolev optical depths (2D array: lines x shells).
        electron_density : ndarray
            Electron densities per shell.
        points : int
            Number of p-ray points.

        Returns
        -------
        L : ndarray
            Integrated luminosities for each frequency.
        intensities_nu_p : ndarray
            Intensities at each p-ray multiplied by p (frequency x p-ray grid).
        """
        return numba_formal_integral(
            self.geometry,
            self.time_explosion,
            self.plasma,
            inner_temperature,
            frequencies,
            att_S_ul,
            Jred_lu,
            Jblue_lu,
            tau_sobolev,
            electron_density,
            points,
        )
