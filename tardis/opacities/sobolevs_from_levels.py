import numpy as np
from numba import njit

from tardis.transport.montecarlo import njit_dict_no_parallel


@njit(**njit_dict_no_parallel)
def calculate_sobolev_opacities_from_level_densities(
    level_density_values,
    line_indices,
    lines_lower_level_index,
    lines_upper_level_index,
    g_lower,
    g_upper,
    meta_stable_upper,
    nlte_lines_mask,
    tau_coefficient,
):
    """Calculate Sobolev optical depths and escape probabilities from level densities.

    Parameters
    ----------
    level_density_values : numpy.ndarray
        Number densities for all atomic levels.
    line_indices : numpy.ndarray, dtype int
        Indices of the lines for which to calculate the Sobolev quantities.
    lines_lower_level_index : numpy.ndarray, dtype int
        Index of the lower level for each line.
    lines_upper_level_index : numpy.ndarray, dtype int
        Index of the upper level for each line.
    g_lower : numpy.ndarray
        Statistical weights of the lower level for each line.
    g_upper : numpy.ndarray
        Statistical weights of the upper level for each line.
    meta_stable_upper : numpy.ndarray, dtype bool
        Mask identifying lines with metastable upper levels.
    nlte_lines_mask : numpy.ndarray, dtype bool
        Mask identifying lines treated with non-local thermodynamic equilibrium
        populations.
    tau_coefficient : numpy.ndarray
        Coefficient multiplying the lower-level density for each line's optical
        depth.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Sobolev optical depths and escape probabilities for the selected lines.
    """
    tau_sobolevs = np.empty(line_indices.shape[0], dtype=np.float64)
    beta_sobolevs = np.empty_like(tau_sobolevs)
    for i in range(line_indices.shape[0]):
        line_index = line_indices[i]
        n_lower = level_density_values[lines_lower_level_index[line_index]]
        n_upper = level_density_values[lines_upper_level_index[line_index]]

        stimulated_emission_factor = 1 - (
            g_lower[line_index] * n_upper / (g_upper[line_index] * n_lower)
        )
        if (
            n_lower == 0.0
            or np.isneginf(stimulated_emission_factor)
            or (
                meta_stable_upper[line_index] and stimulated_emission_factor < 0
            )
            or (nlte_lines_mask[line_index] and stimulated_emission_factor < 0)
        ):
            stimulated_emission_factor = 0.0

        tau_sobolevs[i] = (
            tau_coefficient[line_index] * n_lower * stimulated_emission_factor
        )
        if tau_sobolevs[i] > 1e3:
            beta_sobolevs[i] = tau_sobolevs[i] ** -1
        elif tau_sobolevs[i] < 1e-4:
            beta_sobolevs[i] = 1 - 0.5 * tau_sobolevs[i]
        else:
            beta_sobolevs[i] = (1 - np.exp(-tau_sobolevs[i])) / tau_sobolevs[i]

    return tau_sobolevs, beta_sobolevs
