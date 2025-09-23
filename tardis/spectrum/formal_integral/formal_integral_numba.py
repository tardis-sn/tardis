import numpy as np
from numba import njit, prange
from typing import Tuple
from numpy.typing import NDArray

from tardis.spectrum.formal_integral.base import (
    C_INV,
    BoundsError,
    calculate_impact_parameters,
    intensity_black_body,
)
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.transport.montecarlo import njit_dict, njit_dict_no_parallel
from tardis.transport.montecarlo.configuration.constants import SIGMA_THOMSON


@njit(**njit_dict_no_parallel)
def calculate_intersection_point(
    radius: float, impact_parameter: float, inv_t: float
) -> float:
    """
    Calculate distance to the impact parameter at radius.

    Calculate half of the length of the impact parameter inside a shell
    of radius in terms of unit length (c * t_exp).
    If shell and impact parameter do not intersect, return 0.

    Parameters
    ----------
    radius : float
        Radius of the shell.
    impact_parameter : float
        Distance of the impact parameter to the center of the supernova.
    inv_t : float
        Inverse of the time_explosion, used to normalize to unit length.

    Returns
    -------
    float
        Half the length of the impact parameter inside the shell, or 0 if no intersection.
    """
    if radius > impact_parameter:
        return (
            np.sqrt(radius * radius - impact_parameter * impact_parameter)
            * C_INV
            * inv_t
        )
    else:
        return 0


@njit(**njit_dict_no_parallel)
def populate_intersection_points(
    geometry: NumbaRadial1DGeometry,
    time_explosion: float,
    impact_parameter: float,
    intersection_points: NDArray[np.float64],
    shell_ids: NDArray[np.int64],
) -> int:
    """
    Calculate the intersection points of the impact parameter with each shell.

    Parameters
    ----------
    geometry : object
        Geometry object containing shell radii.
    time_explosion : float
        Time since explosion (seconds).
    impact_parameter : float
        Distance of the integration line to the center.
    intersection_points : ndarray
        Output array to be filled with intersection_point values.
    shell_ids : ndarray
        Output array to be filled with the corresponding shell IDs.

    Returns
    -------
    int
        Number of intersection points found.
    """
    # abbreviations
    r_outer = geometry.r_outer
    N = len(geometry.r_inner)
    inv_t = 1 / time_explosion
    intersection_point = 0
    offset = N

    if impact_parameter <= geometry.r_inner[0]:
        # intersect the photosphere
        for i in range(N):
            intersection_points[i] = 1 - calculate_intersection_point(
                r_outer[i], impact_parameter, inv_t
            )
            shell_ids[i] = i
        return N
    else:
        # no intersection with photosphere
        # that means we intersect each shell twice
        for i in range(N):
            intersection_point = calculate_intersection_point(
                r_outer[i], impact_parameter, inv_t
            )
            if intersection_point == 0:
                continue
            if offset == N:
                offset = i
            # calculate the index in the resulting array
            i_low = N - i - 1  # the far intersection with the shell
            i_up = N + i - 2 * offset  # the nearer intersection with the shell

            intersection_points[i_low] = 1 + intersection_point
            shell_ids[i_low] = i
            intersection_points[i_up] = 1 - intersection_point
            shell_ids[i_up] = i
        return 2 * (N - offset)


@njit(**njit_dict_no_parallel)
def reverse_binary_search(
    x: NDArray[np.float64], x_insert: float, imin: int, imax: int
) -> int:
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
def line_search(
    nu: NDArray[np.float64], nu_insert: float, number_of_lines: int
) -> int:
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


# numba jit
calculate_impact_parameters = njit(
    calculate_impact_parameters, **njit_dict_no_parallel
)
intensity_black_body = njit(intensity_black_body, **njit_dict_no_parallel)


@njit(**njit_dict)
def initialize_formal_integral_inputs(
    frequencies: NDArray[np.float64],
    inner_temperature: float,
    n_impact_parameters: int,
    geometry: NumbaRadial1DGeometry,
    time_explosion: float,
    tau_sobolev: NDArray[np.float64],
) -> Tuple[
    NDArray[np.float64],  # intensities_nu_p
    NDArray[np.float64],  # impact_parameters
    NDArray[np.float64],  # intersection_points
    NDArray[np.int64],  # shell_ids
    NDArray[np.int64],  # n_intersections
    NDArray[np.float64],  # exp_tau_sobolev
]:
    """
    Prepare all arrays and values needed for the loops inside the formal integral.

    Parameters
    ----------
    frequencies : ndarray
        Array of frequency values.
    inner_temperature : float
        Inner boundary temperature.
    n_impact_parameters : int
        Number of impact parameters.
    geometry : object
        Geometry object containing shell radii.
    time_explosion : float
        Time since explosion (seconds).
    tau_sobolev : ndarray
        Sobolev optical depths for each line and shell.

    Returns
    -------
    intensities_nu_p : ndarray
        Intensities at each frequency and impact parameter.
    impact_parameters : ndarray
        Array of impact parameters.
    intersection_points : ndarray
        Array of intersection points for each impact parameter and shell.
    shell_ids : ndarray
        Array of shell IDs for each impact parameter and shell intersection.
    n_intersections : ndarray
        Number of intersections for each impact parameter.
    exp_tau_sobolev : ndarray
        Exponential of negative Sobolev optical depths (flattened).
    """

    n_frequencies = len(frequencies)
    _, size_shell = tau_sobolev.shape
    exp_tau_sobolev = np.exp(-tau_sobolev.T.ravel())
    radius_max = geometry.r_outer[size_shell - 1]
    radius_photosphere = geometry.r_inner[0]

    impact_parameters = np.zeros(n_impact_parameters, dtype=np.float64)
    impact_parameters[::] = calculate_impact_parameters(
        radius_max, n_impact_parameters
    )

    # prepare arrays for returned parameters
    intensities_nu_p = np.zeros(
        (n_frequencies, n_impact_parameters), dtype=np.float64
    )
    intersection_points = np.zeros(
        (n_impact_parameters, 2 * size_shell), dtype=np.float64
    )
    shell_ids = np.zeros((n_impact_parameters, 2 * size_shell), dtype=np.int64)
    n_intersections = np.zeros(n_impact_parameters, dtype=np.int64)

    # loop over frequencies and impact parameters
    for nu_idx in prange(n_frequencies):
        intensities_nu = intensities_nu_p[nu_idx]
        nu = frequencies[nu_idx]
        for impact_parameter_idx in range(1, n_impact_parameters):
            impact_parameter = impact_parameters[impact_parameter_idx]

            # get shell intersections
            n_intersections_p = populate_intersection_points(
                geometry,
                time_explosion,
                impact_parameter,
                intersection_points[impact_parameter_idx],
                shell_ids[impact_parameter_idx],
            )
            n_intersections[impact_parameter_idx] = n_intersections_p
            intersection_point = intersection_points[impact_parameter_idx]

            # if inside the photosphere, set to black body intensity
            # otherwise zero
            if impact_parameter <= radius_photosphere:
                intensities_nu[impact_parameter_idx] = intensity_black_body(
                    nu * intersection_point[0], inner_temperature
                )
            else:
                intensities_nu[impact_parameter_idx] = 0

    return (
        intensities_nu_p,
        impact_parameters,
        intersection_points,
        shell_ids,
        n_intersections,
        exp_tau_sobolev,
    )


@njit(**njit_dict_no_parallel)
def get_electron_scattering_optical_depth(
    escat_optical_depth: float,
    first_contribution_flag: int,
    mean_intensity_red_lu_idx: int,
    intersection_end: float,
    intersection_start: float,
    escat_opacity: float,
    mean_intensity_blue_lu: float,
    mean_intensity_red_lu: float,
    intensities_nu_p: float,
) -> Tuple[float, int, int]:
    """
    Compute the electron scattering optical depth for given segment

    Parameters
    ----------
    escat_optical_depth : float
        Current electron scattering contribution.
    first_contribution_flag : int
        Flag indicating if this is the first contribution (1 if first, 0 otherwise).
    mean_intensity_red_lu_idx : int
        Index for mean_intensity_red_lu.
    intersection_end : float
        Ending intersection point value for the current segment.
    intersection_start : float
        Starting intersection point value for the current segment.
    escat_opacity : float
        Electron scattering opacity.
    mean_intensity_blue_lu : float
        mean intensity of the transition on the blue side for the current segment.
    mean_intensity_red_lu : float
        mean intensity of the transition on the red side for the current segment.
    intensities_nu_p : float
        Intensity at the current frequency and impact parameter.

    Returns
    -------
    escat_optical_depth : float
        Updated electron scattering contribution.
    first_contribution_flag : int
        Updated flag.
    mean_intensity_red_lu_idx : int
        Updated mean_intensity_red_lu index.
    """
    if first_contribution_flag == 1:
        # first contribution to integration
        # NOTE: this treatment of I_nu_b (given
        #   by boundary conditions) is not in Lucy 1999;
        #   should be re-examined carefully
        escat_optical_depth += (
            (intersection_end - intersection_start)
            * escat_opacity
            * (mean_intensity_blue_lu - intensities_nu_p)
        )
        first_contribution_flag = 0
    else:
        # Account for e-scattering, c.f. Eqs 27, 28 in Lucy 1999
        avg_mean_intensity_lu = 0.5 * (
            mean_intensity_red_lu + mean_intensity_blue_lu
        )
        escat_optical_depth += (
            (intersection_end - intersection_start)
            * escat_opacity
            * (avg_mean_intensity_lu - intensities_nu_p)
        )
        # this introduces the necessary offset of one element between
        # the line offset idx
        mean_intensity_red_lu_idx += 1

    return (
        escat_optical_depth,
        first_contribution_flag,
        mean_intensity_red_lu_idx,
    )


@njit(**njit_dict)
def numba_formal_integral(
    geometry: NumbaRadial1DGeometry,
    time_explosion: float,
    plasma,
    inner_temperature: float,
    frequencies: NDArray[np.float64],
    att_S_ul: NDArray[np.float64],
    mean_intensity_red_lu: NDArray[np.float64],
    mean_intensity_blue_lu: NDArray[np.float64],
    tau_sobolev: NDArray[np.float64],
    electron_densities: NDArray[np.float64],
    n_impact_parameters: int,
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
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
    mean_intensity_red_lu : ndarray
        mean intensity of each line transition from upper to lower on the red side for each line and shell.
    mean_intensity_blue_lu : ndarray
        mean intensity of each line transition from upper to lower on the blue side for each line and shell.
    tau_sobolev : ndarray
        Sobolev optical depths for each line and shell.
    electron_densities : ndarray
        Electron densities per shell.
    n_impact_parameters : int
        Number of impact parameters.

    Returns
    -------
    luminosity_densities : ndarray
        Integrated luminosities for each frequency.
    intensities_nu_p : ndarray
        Intensities per frequency and impact parameter
    """
    # Initialize the output which is shared among threads
    n_frequencies = len(frequencies)
    luminosity_densities = np.zeros(n_frequencies, dtype=np.float64)

    radius_max = geometry.r_outer[-1]
    line_list_nu = plasma.line_list_nu
    n_lines = len(line_list_nu)

    # prepare some of the formal integral arrays
    # Inup, p, zs, shell ids, and size of z
    (
        intensities_nu_p,
        impact_parameters,
        intersection_points,
        shell_ids,
        n_intersections,
        exp_tau_sobolev,
    ) = initialize_formal_integral_inputs(
        frequencies,
        inner_temperature,
        n_impact_parameters,
        geometry,
        time_explosion,
        tau_sobolev,
    )

    # loop per frequency
    for nu_idx in prange(n_frequencies):
        intensities_nu = intensities_nu_p[nu_idx]
        nu = frequencies[nu_idx]

        # loop per impact parameter
        for impact_parameter_idx in range(1, n_impact_parameters):
            # get values for this impact parameter
            impact_parameter = impact_parameters[impact_parameter_idx]
            n_intersections_p = n_intersections[impact_parameter_idx]
            intersection_points_p = intersection_points[impact_parameter_idx]
            shell_ids_p = shell_ids[impact_parameter_idx]

            # get the intersection points in frequency space
            # and the corresponding position in the line list
            nu_start = nu * intersection_points_p[0]
            intersection_start = (
                time_explosion / C_INV * (1.0 - intersection_points_p[0])
            )
            idx_nu_start = line_search(line_list_nu, nu_start, n_lines)
            offset = shell_ids_p[0] * n_lines

            nu_ends = nu * intersection_points_p[1:]
            nu_ends_idxs = n_lines - np.searchsorted(
                line_list_nu[::-1], nu_ends, side="right"
            )

            # Initialize "pointers"
            line_idx = int(idx_nu_start)
            line_idx_offset = int(idx_nu_start + offset)
            line_Jred_lu_idx = int(line_idx_offset)

            # flag for first contribution to integration on current impact_parameter
            first_contribution_flag = 1
            escat_optical_depth = 0
            for i in range(n_intersections_p - 1):
                escat_opacity = (
                    electron_densities[int(shell_ids_p[i])] * SIGMA_THOMSON
                )
                nu_end = nu_ends[i]
                nu_end_idx = nu_ends_idxs[i]
                for _ in range(max(nu_end_idx - line_idx, 0)):
                    # calculate e-scattering optical depth to next resonance point
                    intersection_end = (
                        time_explosion
                        / C_INV
                        * (1.0 - line_list_nu[line_idx] / nu)
                    )

                    (
                        escat_optical_depth,
                        first_contribution_flag,
                        line_Jred_lu_idx,
                    ) = get_electron_scattering_optical_depth(
                        escat_optical_depth,
                        first_contribution_flag,
                        line_Jred_lu_idx,
                        intersection_end,
                        intersection_start,
                        escat_opacity,
                        mean_intensity_blue_lu[line_idx_offset],
                        mean_intensity_red_lu[line_Jred_lu_idx],
                        intensities_nu[impact_parameter_idx],
                    )

                    intensities_nu[impact_parameter_idx] += escat_optical_depth
                    # Lucy 1999, Eq 26
                    intensities_nu[impact_parameter_idx] *= exp_tau_sobolev[
                        line_idx_offset
                    ]
                    intensities_nu[impact_parameter_idx] += att_S_ul[
                        line_idx_offset
                    ]

                    # reset e-scattering opacity
                    escat_optical_depth = 0
                    intersection_start = intersection_end

                    line_idx += 1
                    line_idx_offset += 1

                # calculate e-scattering optical depth to grid cell boundary
                avg_mean_intensity_lu = 0.5 * (
                    mean_intensity_red_lu[line_Jred_lu_idx]
                    + mean_intensity_blue_lu[line_idx_offset]
                )
                intersection_end = time_explosion / C_INV * (1.0 - nu_end / nu)
                escat_optical_depth += (
                    (intersection_end - intersection_start)
                    * escat_opacity
                    * (
                        avg_mean_intensity_lu
                        - intensities_nu[impact_parameter_idx]
                    )
                )
                intersection_start = intersection_end

                # advance "pointers" - compute direction on-the-fly
                direction = int((shell_ids_p[i + 1] - shell_ids_p[i]) * n_lines)
                line_idx_offset += direction
                line_Jred_lu_idx += direction
            intensities_nu[impact_parameter_idx] *= impact_parameter
        luminosity_densities[nu_idx] = (
            8
            * np.pi
            * np.pi
            * np.trapezoid(intensities_nu, dx=radius_max / n_impact_parameters)
        )

    return luminosity_densities, intensities_nu_p


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
    n_impact_parameters : int, optional
        Number of impact parameters
    """

    def __init__(
        self,
        geometry: NumbaRadial1DGeometry,
        time_explosion: float,
        plasma,
        n_impact_parameters: int = 1000,
    ):
        self.geometry = geometry
        self.time_explosion = time_explosion
        self.plasma = plasma
        self.n_impact_parameters = n_impact_parameters

    def formal_integral(
        self,
        inner_temperature: float,
        frequencies: NDArray[np.float64],
        att_S_ul: NDArray[np.float64],
        mean_intensity_red_lu: NDArray[np.float64],
        mean_intensity_blue_lu: NDArray[np.float64],
        tau_sobolev: NDArray[np.float64],
        electron_densities: NDArray[np.float64],
        n_impact_parameters: int,
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
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
        mean_intensity_red_lu : ndarray
            mean intensity of each line transition from upper to lower on the red side for each line and shell.
        mean_intensity_blue_lu : ndarray
            mean intensity of each line transition from upper to lower on the blue side for each line and shell.
        tau_sobolev : ndarray
            Sobolev optical depths (2D array: lines x shells).
        electron_densities : ndarray
            Electron densities per shell.
        n_impact_parameters : int
            Number of impact parameters

        Returns
        -------
        luminosity_densities : ndarray
            Integrated luminosities for each frequency.
        intensities_nu_p : ndarray
            Intensities per frequency and impact parameter
        """
        return numba_formal_integral(
            self.geometry,
            self.time_explosion,
            self.plasma,
            inner_temperature,
            frequencies,
            att_S_ul,
            mean_intensity_red_lu,
            mean_intensity_blue_lu,
            tau_sobolev,
            electron_densities,
            n_impact_parameters,
        )
